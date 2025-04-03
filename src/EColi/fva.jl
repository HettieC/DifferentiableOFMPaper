# import AbstractFBCModels.CanonicalModel as CM
# import AbstractFBCModels as A
# import JSONFBCModels as JFBC
# import COBREXA as X
# using JSON, HiGHS
# using ElementaryFluxModes
# using CairoMakie
# import ConstraintTrees as C
# import DifferentiableMetabolism as D

using Distributed
addprocs(14)

@everywhere begin
    import AbstractFBCModels as A
    import AbstractFBCModels.CanonicalModel as CM
    import JSONFBCModels as JFBC
    import COBREXA as X
    using JSON, HiGHS
    import ConstraintTrees as C
    import DifferentiableMetabolism as D
    using ElementaryFluxModes

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
    c = X.enzyme_constrained_flux_balance_constraints(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
    )
end


sol = X.optimized_values(
    c,
    objective=c.objective.value,
    optimizer=HiGHS.Optimizer,
)

objective_bound = X.relative_tolerance_bound(1)
objective = c.objective.value

fva = X.constraints_variability(
    c * :objective_bound^C.Constraint(objective, objective_bound(sol.objective)),
    c.fluxes;
    optimizer=HiGHS.Optimizer,
    workers = workers(),
)

Dict(r => (lb,sol.fluxes[r]) for (r,(lb,ub)) in fva if abs(lb-sol.fluxes[r])>1e-10)
