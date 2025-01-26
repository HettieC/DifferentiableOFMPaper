import AbstractFBCModels.CanonicalModel as CM
import DifferentiableMetabolism as D
import FastDifferentiation as F
const Ex = F.Node
import ConstraintTrees as C
import AbstractFBCModels as A
import JSONFBCModels as JFBC
import COBREXA as X
using JSON, HiGHS
using ElementaryFluxModes
using CairoMakie


flux_zero_tol = 1e-6
gene_zero_tol = 1e-6

model = convert(CM.Model, X.load_model("data/curated_data/EColi/processed_files/iML1515_processed.json"))
model.reactions["EX_glc__D_e"].lower_bound = -1000.0
gene_product_molar_masses = Dict(x => y for (x, y) in JSON.parsefile("data/curated_data/EColi/processed_files/gene_product_molar_mass.json"))
reaction_isozymes = Dict(
    k => Dict(
        "isozyme_$i" => X.Isozyme(d["gene_product_stoichiometry"], d["kcat_forward"], d["kcat_reverse"])
        for (i, d) in enumerate(v)
    )
    for (k, v) in JSON.parsefile("data/curated_data/EColi/processed_files/isozymes.json")
)
reaction_isozymes = Dict(
    rid => Dict("iso" => isozyme[argmax(
        Dict(
            i => iso.kcat_forward /
                 sum([count * gene_product_molar_masses[gene] for (gene, count) in iso.gene_product_stoichiometry])
            for (i, iso) in isozyme
        )
    )])
    for (rid, isozyme) in reaction_isozymes
)

capacity = JSON.parsefile("data/curated_data/EColi/processed_files/capacity.json")
capacity = [(capacity[1][1], string.(capacity[1][2]), capacity[1][3]), (capacity[2][1], string.(capacity[2][2]), capacity[2][3])]

ec_solution = X.enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity,
    optimizer=HiGHS.Optimizer
)

pruned_model, pruned_reaction_isozymes = D.prune_model(
    model,
    ec_solution.fluxes,
    ec_solution.gene_product_amounts,
    reaction_isozymes,
    ec_solution.isozyme_forward_amounts,
    ec_solution.isozyme_reverse_amounts,
    flux_zero_tol,
    gene_zero_tol,
);

parameter_values = Dict(Symbol(x)=>y["iso"].kcat_forward for (x,y) in pruned_reaction_isozymes)
rid_kcat = Dict(k => Ex(Symbol(k)) for (k,_) in parameter_values)
parameter_isozymes = Dict(
    x => Dict(
        "iso" =>  X.IsozymeT{Ex}(
            y["iso"].gene_product_stoichiometry, 
            rid_kcat[Symbol(x)],
            nothing
        )
    )
    for (x,y) in pruned_reaction_isozymes
)

pkm = X.enzyme_constrained_flux_balance_constraints( # kinetic model
    pruned_model;
    reaction_isozymes = parameter_isozymes,
    gene_product_molar_masses,
    capacity,
)
pruned_solution = D.optimized_values(
    pkm,
    parameter_values;
    objective = pkm.objective.value,
    optimizer = HiGHS.Optimizer,
)

# calculate EFMs
N = A.stoichiometry(pruned_model)

# atpm, biomass
atpm_idx = findfirst(x -> x == "ATPM", A.reactions(pruned_model))
biomass_idx = findfirst(x -> x == "BIOMASS_Ec_iML1515_core_75p37M", A.reactions(pruned_model))
fixed_fluxes = [atpm_idx, biomass_idx]
flux_values = [pruned_solution.tree.fluxes["ATPM"], pruned_solution.tree.fluxes["BIOMASS_Ec_iML1515_core_75p37M"]]

OFMs = get_ofms(Matrix(N), fixed_fluxes, flux_values)

OFM_dicts = [
    Dict(A.reactions(pruned_model) .=> OFMs[:,1]),
    Dict(A.reactions(pruned_model) .=> OFMs[:,2])
]
# scale to have one unit flux through biomass
OFM_dicts = [
    Dict(x=>y/OFM_dicts[1]["BIOMASS_Ec_iML1515_core_75p37M"] for (x,y) in OFM_dicts[1]),
    Dict(x=>y/OFM_dicts[2]["BIOMASS_Ec_iML1515_core_75p37M"] for (x,y) in OFM_dicts[2])
]

## calculate lambda
Dict(x => y for (x, y) in OFM_dicts[1] if OFM_dicts[2][x] == 0)
Dict(x => y for (x, y) in OFM_dicts[2] if OFM_dicts[1][x] == 0)
M = [
    OFM_dicts[1]["TALA"] OFM_dicts[2]["TALA"];
    OFM_dicts[1]["ACKr"] OFM_dicts[2]["ACKr"]
]

v = [
    pruned_solution.tree.fluxes["TALA"];
    pruned_solution.tree.fluxes["ACKr"]
]

# rounding causes issues
lambda = inv(M) * v

parameter_values = Dict(Symbol(x)=>first(values(y)).kcat_forward for (x,y) in pruned_reaction_isozymes)
parameters = Ex.(collect(keys(parameter_values)))
rid_pid = Dict(rid => [iso.kcat_forward for (k, iso) in v][1] for (rid, v) in parameter_isozymes)
rid_gcounts = Dict(rid => [v.gene_product_stoichiometry for (k, v) in d][1] for (rid, d) in pruned_reaction_isozymes)

sens_efm = differentiate_efm(OFM_dicts, parameters, rid_pid, parameter_values, rid_gcounts, capacity, gene_product_molar_masses, HiGHS.Optimizer)

param_vals = collect(values(parameter_values))
scaled_sens = Matrix(undef, size(sens_efm, 1), size(sens_efm, 2))
for (i, col) in enumerate(eachcol(sens_efm))
    scaled_sens[:, i] = param_vals[i] .* col ./ lambda
end

order = sortperm(scaled_sens[1, :])
colors = Makie.wong_colors()[5:6]

fig = Figure(; size=(1000, 800), backgroundcolor=:transparent)
data = (
    x=1:length(parameters),
    height=[c[1] > 0 ? c[1] : c[2] for c in eachcol(scaled_sens[:, order])],
    grp=[c[1] > 0 ? 1 : 2 for c in eachcol(scaled_sens[:, order])],
)
ax = Axis(
    fig[1, 1],
    #xlabel = "Parameter",
    ylabel="OFM control coefficient",
    xlabelsize=25,
    xticklabelsize=30,
    ylabelsize=30,
    yticklabelsize=20,
    xgridvisible=false,
    ygridvisible=false,
    xticks=([163, 360], ["Cytosolic\nEnzymes", "Membrane\nEnzymes"]),
    yscale=log10
)
barplot!(
    ax,
    data.x,
    data.height,
    color=colors[data.grp]
)
labels = ["Respirofermentative OFM", "Respiratitive OFM"]
elements = [PolyElement(polycolor=colors[2]), PolyElement(polycolor=colors[1])]
Legend(
    fig[1, 1],
    elements,
    labels,
    position=:ct,
    tellheight=false,
    tellwidth=false,
    halign=:center,
    valign=:top,
    margin=(0, 80, 80, 80),
    labelsize=30,
)
bracket!(ax, 0, 1.6e-6, findlast(x -> x == 2, data.grp), 1.6e-6, style=:curly, orientation=:down)
bracket!(ax, findfirst(x -> x == 1, data.grp), 1.6e-6, length(data.grp), 1.6e-6, style=:curly, orientation=:down)
fig

save("data/plots/ecoli_d_lambda.png", fig)
