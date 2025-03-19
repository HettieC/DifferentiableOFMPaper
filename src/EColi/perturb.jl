import AbstractFBCModels.CanonicalModel as CM
import AbstractFBCModels as A
import JSONFBCModels as JFBC
import COBREXA as X
using JSON, HiGHS
using ElementaryFluxModes
using CairoMakie
import ConstraintTrees as C
import DifferentiableMetabolism as D

flux_zero_tol = 1e-6
gene_zero_tol = 1e-6

# load the model 
model = convert(CM.Model, X.load_model("data/curated_data/EColi/processed_files/iML1515_processed.json"))
model.reactions["EX_glc__D_e"].lower_bound = -1000.0
# set up enzyme constraints
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



# screen through perturbed kcats and see if it makes a non-smooth solution 
res = X.screen(0:0.05:1) do (pgk_kcat)
    c = X.enzyme_constrained_flux_balance_constraints(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
    )
    # increase forward and reverse kcat by 1%
    c.isozyme_flux_forward_balance["PGK"].value.weights[2] = pgk_kcat * reaction_isozymes["PGK"]["iso"].kcat_forward
    c.isozyme_flux_reverse_balance["PGK"].value.weights[3] = pgk_kcat * reaction_isozymes["PGK"]["iso"].kcat_reverse

    sol = X.optimized_values(
        c,
        objective=c.objective.value,
        optimizer=HiGHS.Optimizer,
    )

    # prune the model 
    pruned_model, pruned_reaction_isozymes = D.prune_model(
        model,
        sol.fluxes,
        sol.gene_product_amounts,
        reaction_isozymes,
        sol.isozyme_forward_amounts,
        sol.isozyme_reverse_amounts,
        flux_zero_tol,
        gene_zero_tol,
    )

    # calculate EFMs
    N = A.stoichiometry(pruned_model)

    # atpm, biomass
    atpm_idx = findfirst(x -> x == "ATPM", A.reactions(pruned_model))
    biomass_idx = findfirst(x -> x == "BIOMASS_Ec_iML1515_core_75p37M", A.reactions(pruned_model))
    fixed_fluxes = [atpm_idx, biomass_idx]
    flux_values = [sol.fluxes["ATPM"], sol.fluxes["BIOMASS_Ec_iML1515_core_75p37M"]]

    OFMs = get_ofms(Matrix(N), fixed_fluxes, flux_values)

    OFM_dicts = [
        Dict(A.reactions(pruned_model) .=> OFMs[:, 1]),
        Dict(A.reactions(pruned_model) .=> OFMs[:, 2])
    ]
    # scale to have one unit flux through biomass
    OFM_dicts = [
        Dict(x => y / OFM_dicts[1]["BIOMASS_Ec_iML1515_core_75p37M"] for (x, y) in OFM_dicts[1]),
        Dict(x => y / OFM_dicts[2]["BIOMASS_Ec_iML1515_core_75p37M"] for (x, y) in OFM_dicts[2])
    ]



    return sol
end


f = Figure()
ax = Axis(f[1, 1])
xs = 0:0.05:1
lines!(ax, xs, res)
f





# screen through perturbed kcats and see if OFMs are the same 

rids = collect(keys(reaction_isozymes))[1:100]

res = X.screen(rids) do (rid)
    c = X.enzyme_constrained_flux_balance_constraints(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
    )
    # increase forward and reverse kcat by 1%
    c.isozyme_flux_forward_balance[rid].value.weights[2] = 1.01 * reaction_isozymes["$rid"]["iso"].kcat_forward
    c.isozyme_flux_reverse_balance[rid].value.weights[3] = 1.01 * reaction_isozymes["$rid"]["iso"].kcat_reverse

    sol = X.optimized_values(
        c,
        objective=c.objective.value,
        optimizer=HiGHS.Optimizer,
    )

    # prune the model 
    pruned_model, pruned_reaction_isozymes = D.prune_model(
        model,
        sol.fluxes,
        sol.gene_product_amounts,
        reaction_isozymes,
        sol.isozyme_forward_amounts,
        sol.isozyme_reverse_amounts,
        flux_zero_tol,
        gene_zero_tol,
    )

    # calculate EFMs
    N = A.stoichiometry(pruned_model)

    # atpm, biomass
    atpm_idx = findfirst(x -> x == "ATPM", A.reactions(pruned_model))
    biomass_idx = findfirst(x -> x == "BIOMASS_Ec_iML1515_core_75p37M", A.reactions(pruned_model))
    fixed_fluxes = [atpm_idx, biomass_idx]
    flux_values = [sol.fluxes["ATPM"], sol.fluxes["BIOMASS_Ec_iML1515_core_75p37M"]]

    OFMs = get_ofms(Matrix(N), fixed_fluxes, flux_values)

    OFM_dicts = [
        Dict(A.reactions(pruned_model) .=> OFMs[:, 1]),
        Dict(A.reactions(pruned_model) .=> OFMs[:, 2])
    ]
    # scale to have one unit flux through biomass
    OFM_dicts = [
        Dict(x => y / OFM_dicts[1]["BIOMASS_Ec_iML1515_core_75p37M"] for (x, y) in OFM_dicts[1]),
        Dict(x => y / OFM_dicts[2]["BIOMASS_Ec_iML1515_core_75p37M"] for (x, y) in OFM_dicts[2])
    ]

    return OFM_dicts
end

res2 = [
    [
        Dict(x => y for (x, y) in ofm[1] if y > 1e-5),
        Dict(x => y for (x, y) in ofm[2] if y > 1e-5)
    ]
    for ofm in res
]

# check if any EFMs are different 
diff_ofms = Int64[]
ofm = res2[1]
for (idx, ofms) in enumerate(res2)
    idx == 1 && continue
    if any([x for (x, y) in ofm[1]] .!== [x for (x, y) in ofms[1]])
        push!(diff_ofms, idx)
    end
end


diff_ofms
