using CSV, DataFrames, COBREXA, JSON
import AbstractFBCModels.CanonicalModel as CM
import AbstractFBCModels as A
import MATFBCModels, JSONFBCModels
import AbstractFBCModels: genes, reactions, stoichiometry
include("utils.jl")
organisms = [String(split(x, "_PredictionResults.txt")[1]) for x in readdir("data/data/fungi343species/PredcitedKcat343species") if endswith(x, "_PredictionResults.txt") && !occursin("panmodel",x)]

# make data directories 
!isdir("data/fungi343species") && mkdir("data/fungi343species")
!isdir("data/fungi343species/isozymes") && mkdir("data/fungi343species/isozymes")
!isdir("data/fungi343species/models") && mkdir("data/fungi343species/models")
!isdir("data/fungi343species/molar_mass") && mkdir("data/fungi343species/molar_mass")

# make .mat models into json models, save the isozymes, add fake genes
for organism in organisms
    isfile("data/fungi343species/isozymes/$(organism).json") && isfile("data/fungi343species/models/$organism.json") && continue
    model = convert(CM.Model, load_model("data/data/fungi343species/ssGEMs/$organism.mat"));
    df = select!(
        DataFrame(
            CSV.File("data/data/fungi343species/PredcitedKcat343species/$(organism)_PredictionResults.txt";
                delim='\t',
                types=String)
        ),
        "# rxnID" => "rxnID",
        "Kcat value (substrate first)" => "kcat",
        "genes" => "genes",
    )
    ### collect kcats in correct format, choosing max prediction for each isozyme in each reaction
    ks = Float64[]
    for row in eachrow(df)
        if ismissing(row.kcat) || !contains(row.kcat, r"\d")
            push!(ks, 0.001)
        else
            split_k = split(row.kcat, (',', ';'))
            k = [string(x) for x in split(row.kcat, (',', ';')) if x ∉ ["", "#"]]
            push!(ks, maximum(parse.(Float64, k)))
        end
    end

    kcat_df = DataFrame(
        Reaction=[replace(rxn, "_fwd" => "") for rxn in [rxn[1:6] * replace(rxn[7:end], r"_\d+" => "") for rxn in df.rxnID]],
        Gene=[string.(split(g, ";")) for g in df.genes],
        kcat=ks
    )

    reaction_isozymes = Dict{String,Vector{Isozyme}}()
    for (r, rxn) in model.reactions
        isnothing(rxn.gene_association_dnf) && continue
        reaction_isozymes[r] = Isozyme[]
        for gs in rxn.gene_association_dnf # for each isozyme
            temp_df = filter(row -> row.Reaction == r && all(g -> g ∈ row.Gene, gs), kcat_df)
            kcat_f = isempty(temp_df.kcat) ? 65.0 * 3.6 : temp_df.kcat[1] * 3.6
            if "$(r)_rvs" ∈ kcat_df.Reaction
                temp_df = filter(row -> row.Reaction == r * "_rvs" && all(g -> g ∈ row.Gene, gs), kcat_df)
                kcat_r = isempty(temp_df.kcat) ? 65.0 * 3.6 : temp_df.kcat[1] * 3.6
            else
                kcat_r = 0.001 * 3.6
            end
            if ismissing(kcat_f)
                kcat_f = 65.0 * 3.6
                kcat_r = rxn.lb < 0 ? 65.0 * 3.6 : 0.001 * 3.6
            end
            push!(
                reaction_isozymes[r],
                Isozyme(
                    Dict(
                        g => 1 for g in gs
                    ),
                    kcat_f,
                    kcat_r
                )
            )
        end
    end

    avg_kcat = sum([v[1].kcat_forward for (k, v) in reaction_isozymes if !occursin("transport", A.reaction_name(model, k))]) /
               length([v[1].kcat_forward for (k, v) in reaction_isozymes if !occursin("transport", A.reaction_name(model, k))])

    avg_transport = sum([v[1].kcat_forward for (k, v) in reaction_isozymes if occursin("transport", A.reaction_name(model, k))]) /
                    length([v[1].kcat_forward for (k, v) in reaction_isozymes if occursin("transport", A.reaction_name(model, k))])

    # add fake isozymes to reactions with no gene 
    forward_isozyme = Isozyme(
        Dict("s0001" => 1),
        avg_kcat,
        0.01 * 3.6
    )
    reversible_isozyme = Isozyme(
        Dict("s0001" => 1),
        avg_kcat,
        avg_kcat
    )
    fake_transport = Isozyme(
        Dict("s0002" => 1),
        avg_transport,
        avg_transport
    )

    model.genes["s0001"] = CM.Gene(; name="s0001")
    model.genes["s0002"] = CM.Gene(; name="s0002")
    for (rid, rxn) in model.reactions
        (haskey(reaction_isozymes, rid) || rid == "r_4041" || rid == "r_4046" || occursin("exchange", rxn.name)) && continue
        if contains(rxn.name, "transport")
            reaction_isozymes[rid] = [fake_transport]
            model.reactions[rid].gene_association_dnf = [["s0002"]]
        elseif model.reactions[rid].lower_bound < 0
            reaction_isozymes[rid] = [reversible_isozyme]
            model.reactions[rid].gene_association_dnf = [["s0001"]]
        else
            reaction_isozymes[rid] = [forward_isozyme]
            model.reactions[rid].gene_association_dnf = [["s0001"]]
        end
    end

    open("data/fungi343species/isozymes/$(organism).json", "w") do io
        JSON.print(io, reaction_isozymes)
    end
    save_model(convert(JSONFBCModels.JSONFBCModel, model), "data/fungi343species/models/$organism.json");
end

## save the gene product molar mass and mass group
for organism in organisms
    model = convert(CM.Model, load_model("data/fungi343species/models/$organism.json"));
    gpmm = get_gene_product_molar_mass(organism, model)
    avg_mw = sum(values(gpmm)) / length(gpmm)
    for x in [x for (x, y) in model.genes if !haskey(gpmm, x)]
        gpmm[x] = avg_mw
    end
    open("data/fungi343species/molar_mass/$organism.json", "w") do io
        JSON.print(io, gpmm)
    end
end
