using CairoMakie, CSV, MAT, JSON, ColorSchemes
using DataFrames, COBREXA, Statistics, MATFBCModels
import AbstractFBCModels.CanonicalModel as CM


### Plot how sensitive biomass is to reactions in different subsystems
# load pan model to get the subsystems

m = matread("data/fungi343species/ssGEMs/panmodel.mat")["panmodel"]
model = convert(CM.Model, load_model("data/fungi343species/ssGEMs/panmodel.mat"));
reaction_subsystems = Dict(m["rxns"] .=> vec.(m["subSystems"]))
reaction_subsystems = Dict(x => join([z for z in y], "; ") for (x, y) in reaction_subsystems)
pathways = Dict(x => y for (x, y) in JSON.parsefile("data/fungi343species/database/kegg_pathways.json"))

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

organisms = [String(split(x, ".csv")[1]) for x in readdir("data/fungi343species/fixed_EC_sensitivities")]

# for each model, save the avg sens of biomass to reactions in each subsystem 
subsystem_sens = Dict(p => Float64[] for (p, v) in pathways)
for organism in organisms
    bio_df = filter!(
        row -> row.vid == "(:fluxes, :r_2111)",
        DataFrame(CSV.File("data/fungi343species/fixed_EC_sensitivities/$organism.csv")),
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


save("data/plots/yeasts_subs.png", f)

# for each model, save the percentage of reactions per subsystem in the core model
all_rxns = [string.(vec(matread("data/fungi343species/ssGEMs/$organism.mat")["reducedModel"]["rxns"])) for organism in organisms]
pathway_rxns = Dict(p => String[] for (p, q) in pathways)
pathway_rxns["Other"] = String[]
for (p, q) in pathway_rxns
    for (k, v) in simplified_rxn_subsystems
        if p ∈ v
            push!(pathway_rxns[p], k)
        end
    end
end


core_rxns = String[]
model_rxns = Dict{String,Vector{String}}()
for organism in organisms
    model_rxns[organism] = sort(vec(matread("data/fungi343species/ssGEMs/$organism.mat")["reducedModel"]["rxns"]))
end

# model with fewest rxns is Schizosaccharomyces_pombe
core_rxns = model_rxns["Schizosaccharomyces_pombe"]
delete!(model_rxns, "Schizosaccharomyces_pombe")
for (x, y) in model_rxns
    core_rxns = intersect(y, core_rxns)
end

core_pathways = Dict(x => filter(z -> z ∈ core_rxns, y) for (x, y) in pathway_rxns)

model_subsystem_conservation = Dict(model => Dict(sub => length(core_pathways[sub]) / length([r for r in rxns if r ∈ ps]) for (sub, ps) in pway_rxn) for (model, rxns) in model_rxns)



conservation = Dict(p => length(v) / length(pathway_rxns[p]) for (p, v) in core_pathways)
open("data/fungi343species/conservation.json", "w") do io
    JSON.print(io, conservation)
end
open("data/fungi343species/conserved_rxns.json", "w") do io
    JSON.print(io, core_pathways)
end
cons_vals = [y for (x, y) in conservation]


# conservation per model 
model_pathways = Dict(
    k => Dict(x => [y for y in pway_rxn[x] if y ∈ v] for (x, y) in pway_rxn if x != "Other")
    for (k, v) in model_rxns
)

# TODO finish!
model_pathway_conservation = Dict(
    k => Dict(x => length(y) / length(v[x]) for (x, y) in core_pathways) for (k, v) in model_pathways
)


############## colour the plot by how many of the reactions in the model are in the core model
sens_keys = [replace(x, (" " => "\n")) for (x, y) in subsystem_sens]
sens_vals = [y for (x, y) in subsystem_sens]
order = sortperm(sens_vals, by=x -> median(x))

cats = vcat([repeat([i], length(organisms)) for i in 1:10]...)
vals = log.(vcat(sens_vals[order]...) .+ 1e-10)
xlabs = sens_keys[order]

conservation = JSON.parsefile("data/fungi343species/conservation.json")
labels = [x for (x, y) in subsystem_sens][order]
cons_vals = [conservation[x] for x in labels]
colors = reverse(colorschemes[:plasma])
norm_cons = [(x - minimum(cons_vals)) / (maximum(cons_vals) - minimum(cons_vals)) for x in cons_vals]
eyes = Int.(floor.([255 * x + 1 for x in norm_cons]))
cols = colors[vcat([repeat([x], 343) for x in eyes]...)]


# make the plot
f = Figure(; size=(1600, 1000), background_color=:transparent,);
ax = Axis(
    f[1, 1],
    xlabel="Subsystem",
    ylabel="Mean sensitivity of biomass to subsystem enzymes",
    xlabelsize=35,
    ylabelsize=35,
    titlesize=30,
    xticklabelsize=30,
    yticklabelsize=30,
    backgroundcolor=:transparent,
    xticks=(1:length(subsystem_sens), xlabs),
    ytickformat=vals -> [rich("10", superscript("$(trunc(Int,v))")) for v in vals],
    xticklabelrotation=0.8pi / 3
)
rainclouds!(
    ax,
    cats,
    vals;
    #cloud_width=1.8,
    color=cols,
    #gap=0.05,
    plot_boxplots=false,
    show_median=false
)
ylims!(-13, 0)
Colorbar(
    f[1, 2],
    limits=(minimum(values(conservation)), 0.6),
    colormap=reverse(colorschemes[:plasma]),
    size=25,
    label="Proportion of total reactions conserved among models",
    labelsize=30,
    ticklabelsize=25,
    ticks=0.2:0.05:0.65
)
f

