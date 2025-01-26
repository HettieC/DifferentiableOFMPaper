import AbstractFBCModels.CanonicalModel as CM
import DifferentiableMetabolism as D
import FastDifferentiation as F
const Ex = F.Node
import ConstraintTrees as C
import AbstractFBCModels as A
import COBREXA as X
using JSON, HiGHS
using JSONFBCModels 

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

for (r, rxn) in model.reactions
    if rxn.lower_bound < 0 && ec_solution.fluxes[r] < 0
        println(r)
    end
end

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

parameters = collect(keys(parameter_values))

pkm_kkt, vids = D.differentiate_prepare_kkt(pkm, pkm.objective.value, parameters)

sens = D.differentiate_solution(
    pkm_kkt,
    pruned_solution.primal_values,
    pruned_solution.equality_dual_values,
    pruned_solution.inequality_dual_values,
    parameter_values,
    scale = true, # unitless sensitivities
)
