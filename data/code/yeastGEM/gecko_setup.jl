using COBREXA, JSON, DifferentiableMetabolism
import AbstractFBCModels.CanonicalModel as CM
import AbstractFBCModels as A
using CSV, DataFrames, JSONFBCModels
# load the chalmers yeast9 model
model = convert(CM.Model,load_model("data/curated_data/yeastGEM/yeastGEMcurated.json"))

## kcat data
df = convert.(String, DataFrame(CSV.File("data/curated_data/yeastGEM/reaction_data.csv")))
kcat_df = DataFrame(CSV.File("data/curated_data/yeastGEM/turnup_kcats.csv"))[:, 2:end]
big_df = unique!(innerjoin(kcat_df, df, on=[:enzyme, :substrates, :products]))[:, [8, 9, 7]]

protein_stoichiometry = JSON.parsefile("data/curated_data/yeastGEM/protein_stoichiometry.json")

##### make isozymes for all reactions 
reaction_isozymes = Dict{String,Vector{Isozyme}}()
not_found = []
for (r, rxn) in model.reactions
    isnothing(rxn.gene_association_dnf) && continue
    reaction_isozymes[r] = Isozyme[]
    for (i,gs) in enumerate(rxn.gene_association_dnf) # for each isozyme
        temp_gs = [split(g, '_')[1] for g in gs]
        if all(g -> g ∉ big_df.Gene, temp_gs)
            #= if none of the genes in an isozyme are in the dataframe
            then use ~65.0 for the kcats =#
            push!(
                reaction_isozymes[r],
                Isozyme(
                    Dict(
                        k => v for (k,v) in zip(
                            gs,protein_stoichiometry[r][i]
                        )    
                    ),
                    65.0 * 3.6, # need 1/h units
                    rxn.lower_bound < 0 ? 65.0 * 3.6 : 0.01 * 3.6
                )
            )
        else
            new_df = filter(row -> row.Reaction == r && row.Gene ∈ temp_gs, big_df)
            if size(new_df,1) == 0 
                push!(not_found,[r,gs])
            else
                kcat_f = maximum(filter(row -> row.Reaction == r && row.Gene ∈ temp_gs, big_df)[:, 3]) * 3.6
                # otherwise take the maximum kcat for the isozyme, "gs"
                if "$(r)_rev" ∈ big_df.Reaction
                    kcat_r = maximum(filter(row -> row.Reaction == "$(r)_rev" && row.Gene ∈ temp_gs, big_df)[:, 3]) * 3.6
                else
                    kcat_r = 0.01 * 3.6
                end
                if ismissing(kcat_f)
                    kcat_f = 65.0 * 3.6
                    kcat_r = rxn.lower_bound < 0 ? 65.0 * 3.6 : 0.01 * 3.6
                end
                push!(
                    reaction_isozymes[r],
                    Isozyme(
                        Dict(
                            k => v for (k,v) in zip(
                                gs,protein_stoichiometry[r][i]
                            )    
                        ),
                        kcat_f,
                        kcat_r
                    )
                )
            end
        end
    end
end


# add fake isozymes to reactions with no gene 
forward_isozyme = Isozyme(
        Dict("s0001"=>1),
        45.0 * 3.6,
        0.01 * 3.6
)
reversible_isozyme = Isozyme(
    Dict("s0001"=>1),
    45.0 * 3.6,
    45.0 * 3.6
)
fake_transport = Isozyme(
    Dict("s0002"=>1),
    35.0 * 3.6,
    35.0 * 3.6
)
model.genes["s0001"] = CM.Gene()
model.genes["s0002"] = CM.Gene()
for rid in [r for r in A.reactions(model) if !haskey(reaction_isozymes, r) && !contains(A.reaction_name(model, r), "exchange") && r != "r_4041" && r != "r_4046"]
    if !isnothing(A.reaction_gene_association_dnf(model, rid))
        reaction_isozymes[rid] = [
            Isozyme(
                Dict(
                    k => v for (k, v) in zip(
                        reaction_gene_association(model, rid)[i],
                        haskey(reaction_protein_stoichiometry, rid) ? reaction_protein_stoichiometry[rid][i] : 1.0,
                    )
                    if !isnothing(reaction_gene_association(model, rid))
                ),
                contains(A.reaction_name(model, rid), "transport") ? 45.0 * 3.6 : 65.0 * 3.6,
                model.reactions[rid].lb < 0 ? 45.0 * 3.6 : 0.01 * 3.6
            ) for i in 1:length(reaction_gene_association(model, rid))
        ]
    elseif contains(A.reaction_name(model, rid), "transport")
        reaction_isozymes[rid] = [fake_transport]
        global model.reactions[rid].gene_association_dnf = [["s0002"]]
    elseif model.reactions[rid].lower_bound < 0 
        reaction_isozymes[rid] = [reversible_isozyme]
        global model.reactions[rid].gene_association_dnf = [["s0001"]]
    else
        reaction_isozymes[rid] = [forward_isozyme]
        global model.reactions[rid].gene_association_dnf = [["s0001"]]
    end
end
open("data/curated_data/yeastGEM/reaction_isozymes.json","w") do io 
    JSON.print(io, reaction_isozymes)
end

## get gene molar masses
AA_df = DataFrame(CSV.File("data/data/yeastGEM/databases/swissprot.tsv"))
gene_product_molar_mass = Dict{String,Float64}()
for row in eachrow(AA_df)
    gs = split(row.gene_id)
    for g in gs
        if g ∈ A.genes(model)
            gene_product_molar_mass[g] = row.MW / 1000.0 # in Da (g/mol) originally, make kg/mol (kDa)
        end
    end
end
avg_mw = sum(values(gene_product_molar_mass)) / length(gene_product_molar_mass)
for x in [x for x in A.genes(model) if !haskey(gene_product_molar_mass, x)]
    gene_product_molar_mass[x] = avg_mw + rand([-1, 1]) * rand()
end
open("data/curated_data/yeastGEM/gene_product_molar_mass.json","w") do io
    JSON.print(io,gene_product_molar_mass)
end

##### Get proteomics data
### total protein g per gram dry weigh biomass: 0.396, doi:https://doi.org/10.1021/jf0400821
### 13% of yeast proteins are mitochondrial http://www.genesdev.org/cgi/doi/10.1101/gad.970902.
### ∼17% of the ∼6000 different proteins synthesized in a yeast cell.  https://doi.org/10.1091%2Fmbc.E05-08-0740
### 44.9% of mitochondrial proteins are metabolic, mitochondrial proteome https://doi.org/10.1016/j.celrep.2017.06.014
### file mcp.M115.054288-4.xlsx taken from https://doi.org/10.1074/mcp.M115.054288
### unified proteomics dataset https://doi.org/10.1016/j.cels.2017.12.004
### dry weight of 1 cell yeast https://doi.org/10.1371/journal.pone.0188388
proteomics_df = DataFrame(CSV.File("data/data/yeastGEM/proteomics/1-s2.0-S240547121730546X-mmc4.csv"))
proteomics_quant = Dict(Pair.(proteomics_df.SystematicName,proteomics_df.MedianMoleculesPerCell))
avg_abundance = floor(sum([y for (x,y) in proteomics_quant if !ismissing(y)])/length([x for (x,y) in proteomics_quant if !ismissing(y)]))
for (k,v) in proteomics_quant
    if ismissing(v)
        proteomics_quant[k] = floor(avg_abundance/100) # give a very small abundance to proteins with missing values
    end
end

N = 6.02214076E23 # Avogadro constant 
# grams of each protein per gDW cells
prot_g_gDW = Dict(
    x => n * (1/N) * gene_product_molar_mass[x] * 1000 * 4.765E11 
    for (x,n) in proteomics_quant if haskey(gene_product_molar_mass,x)
)
open("data/curated_data/yeastGEM/prot_g_gDW.json","w") do io 
    JSON.print(io, prot_g_gDW)
end

mitochondria_df = DataFrame(CSV.File("data/data/yeastGEM/proteomics/1-s2.0-S2211124717308112-mmc4.csv")) # mitochondrial proteome

mitochondria = String[]
for (k,v) in model.reactions 
    isnothing(v.gene_association_dnf) && continue 
    if any(((m,s),) -> A.metabolite_compartment(model,m) == "m",v.stoichiometry)
        append!(mitochondria,vcat(v.gene_association_dnf...))
    end
end
unique!(mitochondria)


# NOTE: some reactions have gene_association_dnfs where one is mitochondrial one isn't, in this case say that both are mitochondrial
gene_product_mass_group = Dict{String,String}()
for (k,v) in model.reactions 
    isnothing(v.gene_association_dnf) && continue
    if any(g -> g ∈ mitochondria, vcat(v.gene_association_dnf...))
        for g in vcat(v.gene_association_dnf...)
            gene_product_mass_group[g] = "mitochondria"
        end 
    else 
        for g in vcat(v.gene_association_dnf...)
            haskey(gene_product_mass_group,g) && continue
            gene_product_mass_group[g] = "other"
        end
    end
end
gene_product_mass_group["s0001"] = "other"
gene_product_mass_group["s0002"] = "other"


# use the consolidated dataset to get the abundance of mitochondrial proteins, so that it's more comparable to the other cell proteins
mito_total = sum(prot_g_gDW[k] for (k,v) in gene_product_mass_group if v == "mitochondria" && haskey(prot_g_gDW,k))
other_total = sum(prot_g_gDW[k] for (k,v) in gene_product_mass_group if v == "other" && haskey(prot_g_gDW,k))

gene_product_mass_group_bound = Dict("mitochondria"=>mito_total*1000,"other"=>other_total*1000)

capacity = [
    (
        "mitochondria",
        [x for (x,y) in gene_product_mass_group if y == "mitochondria"],
        gene_product_mass_group_bound["mitochondria"]
    ),
    (
        "other",
        [x for (x,y) in gene_product_mass_group if y == "other"],
        gene_product_mass_group_bound["other"]
    )
]

open("data/curated_data/yeastGEM/capacity.json","w") do io 
    JSON.print(io,capacity)
end
