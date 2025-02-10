using COBREXA, JSONFBCModels
import AbstractFBCModels.CanonicalModel as CM
import AbstractFBCModels: genes, reactions, stoichiometry
import AbstractFBCModels as A
using JSON, CSV, DataFrames
include("utils.jl")
model = convert(CM.Model,load_model("data/EColi/iML1515.json"))

pid_bid = Dict()
open(joinpath("data", "EColi", "e_coli.tab")) do io
    firstline = true
    for ln in eachline(io)
        firstline && (firstline = false; continue)
        prts = split(ln, "\t")
        prts[3] == "" && continue
        pid = prts[1]
        bnum = first(split(prts[3], " "))
        pid_bid[pid] = bnum
    end
end

# load protein data
proteome_data = Dict()
subcellular_location = Dict{String,Vector{String}}()
open("data/EColi/uniprotkb_e_coli.tsv") do io
    firstline = true
    for ln in eachline(io)
        firstline && (firstline = false; continue)
        prts = split(ln, "\t")
        !haskey(pid_bid,prts[1]) && continue
    
        gene = pid_bid[prts[1]]
        mass = parse(Float64, replace(prts[2], "," => "")) / 1000.0 # in Da originally
        if prts[4] == ""
            subunit = ""
        else
            subunit =
                split(replace(replace(prts[4], "SUBUNIT: " => ""), "." => ""), " ")
            subunit = subunit[1]
            if !contains(subunit, "mer")
                subunit = ""
            end
        end
        proteome_data[gene] = (mass, subunit)
        if prts[5] == ""
            subcellular_location[gene] = ["uncategorised"]
        else
            subs = string.(split(split(prts[5],"LOCATION: ")[2],"; "))
            subcellular_location[gene] = rstrip.(first.(split.(subs," {")),'.')
        end
    end
end
membrane_proteins = [x for (x,y) in subcellular_location if x ∈ genes(model) && any(z -> occursin("membrane",lowercase(z)),y) && "Cytoplasm" ∉ y]


# get complex data
complex_data = Dict()
open(joinpath("data", "EColi", "83333.tsv")) do io
    firstline = true
    for ln in eachline(io)
        firstline && (firstline = false; continue)
        prts = split(ln, "\t")
        id = prts[1]
        raw_stoich = prts[5]
        stoich = Dict{String,Float64}()
        for elem in split(raw_stoich, "|")
            if startswith(elem, "P") || startswith(elem, "Q")
                pid_stoich = split(elem, "(")
                pid = first(split(first(pid_stoich), "-"))
                stoich_num = first(split(last(pid_stoich), ")"))
                !haskey(pid_bid, pid) && continue
                bid = pid_bid[pid]
                stoich[bid] = parse(Float64, stoich_num)
            end
        end
        isempty(stoich) && continue
        complex_data[id] = stoich
    end
end


#! remove really low abundance metabolites from biomass - numerical precision issues
for (k, v) in model.reactions["BIOMASS_Ec_iML1515_core_75p37M"].stoichiometry
    if abs(v) < 1e-4 # cutoff
        println("Deleted ", k, " with coefficient ", v)
        delete!(model.reactions["BIOMASS_Ec_iML1515_core_75p37M"].stoichiometry, k)
    end
end
#! delete other biomass reaction
delete!(model.reactions, "BIOMASS_Ec_iML1515_WT_75p37M")

reaction_kcats = JSON.parsefile("data/EColi/e_coli_kcats.json")
gene_product_molar_mass = Dict(
    k => v[1] for (k,v) in proteome_data
)
avg_mw = sum(values(gene_product_molar_mass)) / length(gene_product_molar_mass);
gene_product_molar_mass["s0001"] = avg_mw;
gene_product_molar_mass["s0002"] = avg_mw;
gene_product_molar_mass["b1692"] = avg_mw;
open("data/curated_data/EColi/processed_files/gene_product_molar_mass.json","w") do io 
    JSON.print(io,gene_product_molar_mass)
end
proteome_data["s0001"] = (avg_mw, "")
proteome_data["s0002"] = (avg_mw, "")
reaction_protein_stoichiometry = get_protein_stoichiometry!(model,complex_data, proteome_data)


# create isozymes
reaction_isozymes = Dict{String,Vector{Isozyme}}()
for rid in [r for r in keys(reaction_kcats) if r ∈ reactions(model)]
    isnothing(A.reaction_gene_association_dnf(model,rid)) && continue
    reaction_isozymes[rid] = [
        Isozyme(
            Dict(
                k => v for (k,v) in zip(
                    A.reaction_gene_association_dnf(model,rid)[i],
                    reaction_protein_stoichiometry[rid][i],
                )
            ),
            reaction_kcats[rid]*3.6,
            reaction_kcats[rid]*3.6,
        ) for i = 1:length(reaction_kcats[rid])
    ]
end

avg_kcat = sum([v[1].kcat_forward for (k,v) in reaction_isozymes if !occursin("transport",A.reaction_name(model,k))]) /
    length([v[1].kcat_forward for (k,v) in reaction_isozymes if !occursin("transport",A.reaction_name(model,k))])

avg_transport = sum([v[1].kcat_forward for (k,v) in reaction_isozymes if occursin("transport",A.reaction_name(model,k))]) /
    length([v[1].kcat_forward for (k,v) in reaction_isozymes if occursin("transport",A.reaction_name(model,k))])

# add fake isozymes to reactions with no gene 
forward_isozyme = Isozyme(
        Dict("s0001"=>1),
        avg_kcat,
        0.01 * 3.6
)
reversible_isozyme = Isozyme(
    Dict("s0001"=>1),
    avg_kcat,
    avg_kcat
)
fake_transport = Isozyme(
    Dict("s0002"=>1),
    avg_transport,
    avg_transport
)
model.genes["s0001"] = CM.Gene(;name="s0001")
model.genes["s0002"] = CM.Gene(;name="s0002")

# add fake isozymes to reactions with no gene 
for rid in [r for r in reactions(model) if !haskey(reaction_isozymes, r) && !contains(A.reaction_name(model, r), "exchange") && r != "BIOMASS_Ec_iML1515_core_75p37M" && r != "BIOMASS_Ec_iML1515_WT_75p37M"]
    if !isnothing(A.reaction_gene_association_dnf(model, rid))
        reaction_isozymes[rid] = [
            Isozyme(
                Dict(
                    k => v for (k, v) in zip(
                        A.reaction_gene_association_dnf(model, rid)[i],
                        haskey(reaction_protein_stoichiometry, rid) ? reaction_protein_stoichiometry[rid][i] : 1.0,
                    )
                    if !isnothing(A.reaction_gene_association_dnf(model, rid))
                ),
                contains(A.reaction_name(model, rid), "transport") ? 45.0 * 3.6 : 65.0 * 3.6,
                model.reactions[rid].lower_bound < 0 ? 45.0 * 3.6 : 0.01 * 3.6
            ) for i in 1:length(A.reaction_gene_association_dnf(model, rid))
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

open("data/curated_data/EColi/processed_files/isozymes.json","w") do io 
    JSON.print(io,reaction_isozymes)
end

# create enzyme capacity constraints for cytosol and periplasm
gene_product_bounds(gid) = (0, 1000000.0)
gene_product_mass_group = Dict(g => "Periplasm" ∈ subcellular_location[g] ? "periplasm" : "other" for g in genes(model) if haskey(subcellular_location,g))
gene_product_mass_group = Dict(x => x ∈ membrane_proteins ? "membrane" : "other" for x in genes(model))
gene_product_mass_group["s0001"] = "other"
gene_product_mass_group["s0002"] = "other"
gene_product_mass_group["b1692"] = "other"


save_model(convert(JSONFBCModels.JSONFBCModel,model),"data/curated_data/EColi/processed_files/iML1515_processed.json")

mass_percentage = DataFrame(CSV.File("data/EColi/41587_2016_BFnbt3418_MOESM18_ESM_percentage_mass_S14.csv"))

filter!(row->!occursin("Standard",row.Subcellular_location) && !occursin("Ratio",row.Subcellular_location),mass_percentage)

gene_product_mass_group_bound = Dict(
    "membrane" => 5*sum([row.Glucose for row in eachrow(mass_percentage) if occursin("embrane",row.Subcellular_location) || row.Subcellular_location == "Periplasm"]),
    "other" => 5*sum([row.Glucose for row in eachrow(mass_percentage) if row.Subcellular_location == "Cytoplasm" || row.Subcellular_location == "Secreted"])
)

capacity = [
    (
        "membrane",
        [x for (x,y) in gene_product_mass_group if y=="membrane"],
        gene_product_mass_group_bound["membrane"]
    ),
    (
        "other",
        [x for (x,y) in gene_product_mass_group if y=="other"],
        gene_product_mass_group_bound["other"]
    )
]
open("data/curated_data/EColi/processed_files/capacity.json","w") do io 
    JSON.print(io,capacity)
end
