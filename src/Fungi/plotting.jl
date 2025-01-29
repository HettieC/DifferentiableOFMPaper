using CairoMakie, CSV, MAT, JSON, ColorSchemes
using DataFrames, COBREXA, Statistics, MATFBCModels
import AbstractFBCModels.CanonicalModel as CM

### Plot how sensitive biomass is to reactions in different subsystems
# load pan model to get the subsystems

m = matread("data/data/fungi343species/ssGEMs/panmodel.mat")["panmodel"]
model = convert(CM.Model, load_model("data/data/fungi343species/ssGEMs/panmodel.mat"));
reaction_subsystems = Dict(m["rxns"] .=> vec.(m["subSystems"]))
reaction_subsystems = Dict(x => join([z for z in y], "; ") for (x, y) in reaction_subsystems)
pathways = Dict(x => y for (x, y) in JSON.parsefile("data/curated_data/kegg_pathways.json"))

##### keep reactions having more than one subsystem
simplified_rxn_subsystems = Dict(r => String[] for (r, v) in reaction_subsystems)
for (r, subs) in reaction_subsystems
    if any([x ∈ pathways["Carbohydrate"] for x in split(subs, "; ")])
        push!(simplified_rxn_subsystems[r], "Carbohydrate")
    end
    if any([x ∈ pathways["Energy"] for x in split(subs, "; ")])
        push!(simplified_rxn_subsystems[r], "Energy")
    end
    if any([x ∈ pathways["Lipid"] for x in split(subs, "; ")])
        push!(simplified_rxn_subsystems[r], "Lipid")
    end
    if any([x ∈ pathways["Nucleotide"] for x in split(subs, "; ")])
        push!(simplified_rxn_subsystems[r], "Nucleotide")
    end
    if any([x ∈ pathways["Amino acid"] for x in split(subs, "; ")])
        push!(simplified_rxn_subsystems[r], "Amino acid")
    end
    if any([x ∈ pathways["Cofactors, vitamins"] for x in split(subs, "; ")])
        push!(simplified_rxn_subsystems[r], "Cofactors, vitamins")
    end
    if any([x ∈ pathways["Glycan"] for x in split(subs, "; ")])
        push!(simplified_rxn_subsystems[r], "Glycan")
    end
    if any([x ∈ pathways["Terpenoids, polyketides"] for x in split(subs, "; ")])
        push!(simplified_rxn_subsystems[r], "Terpenoids, polyketides")
    end
    if any([x ∈ pathways["Genetic"] for x in split(subs, "; ")])
        push!(simplified_rxn_subsystems[r], "Genetic")
    end
    if any([x ∈ pathways["Cellular processes"] for x in split(subs, "; ")])
        push!(simplified_rxn_subsystems[r], "Cellular processes")
    end
    if isempty(simplified_rxn_subsystems[r])
        simplified_rxn_subsystems[r] = ["Other"]
    end
end

pway_rxn = Dict(p => [x for (x, y) in simplified_rxn_subsystems if p ∈ y] for p in collect(keys(pathways)))

organisms = [String(split(x, ".csv")[1]) for x in readdir("data/data/fungi343species/fixed_EC_sensitivities")]

# for each model, save the avg sens of biomass to reactions in each subsystem 
subsystem_sens = Dict(p => Float64[] for (p, v) in pathways)
for organism in organisms
    bio_df = filter!(
        row -> row.vid == "(:fluxes, :r_2111)",
        DataFrame(CSV.File("data/data/fungi343species/fixed_EC_sensitivities/$organism.csv")),
    )
    bio_sens = Dict(names(bio_df)[2:end] .=> Vector(bio_df[1, 2:end]))

    for (p, vals) in subsystem_sens
        push!(
            vals,
            mean(
                v for (k, v) in bio_sens
                if haskey(simplified_rxn_subsystems, k)[1]
                &&
                p ∈ simplified_rxn_subsystems[k]
            )
        )
    end
end

sens_keys = [replace(x, (" " => "\n")) for (x, y) in subsystem_sens]
sens_vals = [y for (x, y) in subsystem_sens]
order = sortperm(sens_vals, by=x -> median(x))

cats = vcat([repeat([i], length(organisms)) for i in 1:10]...)
vals = log.(vcat(sens_vals[order]...) .+ 1e-10)
xlabs = sens_keys[order]

colors = (ColorSchemes.rainbow)[[1, 2, 4, 6, 8, 10, 12, 14, 15, 17]]
colors = vcat([repeat([color], length(organisms)) for color in colors]...)
# make the plot
f = Figure(; size=(1000, 700), backgroundcolor=:transparent);
ax = Axis(
    f[1, 1],
    backgroundcolor=:transparent,
    xlabel="Subsystem",
    ylabel="Mean sensitivity of biomass to subsystem enzymes",
    xlabelsize=35,
    ylabelsize=35,
    titlesize=30,
    xticklabelsize=25,
    yticklabelsize=25,
    xticks=(1:length(subsystem_sens), xlabs),
    ytickformat=vals -> [rich("10", superscript("$(trunc(Int,v))")) for v in vals],
    xticklabelrotation=pi / 3,
    ygridvisible=false,
)
rainclouds!(
    ax,
    cats,
    vals;
    cloud_width=1.8,
    gap=0.05,
    plot_boxplots=false,
    show_median=true,
    color=colors
)
ylims!(-13, 0)
f
