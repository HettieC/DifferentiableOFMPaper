using Distributed
addprocs(14)

# # e coli fva
# @everywhere begin
#     import AbstractFBCModels as A
#     import AbstractFBCModels.CanonicalModel as CM
#     import JSONFBCModels as JFBC
#     import COBREXA as X
#     using JSON, HiGHS
#     import ConstraintTrees as C
#     import DifferentiableMetabolism as D
#     using ElementaryFluxModes

#     flux_zero_tol = 1e-6
#     gene_zero_tol = 1e-6

#     # load the model 
#     model = convert(CM.Model, X.load_model("data/curated_data/EColi/processed_files/iML1515_processed.json"))
#     model.reactions["EX_glc__D_e"].lower_bound = -1000.0
#     # set up enzyme constraints
#     gene_product_molar_masses = Dict(x => y for (x, y) in JSON.parsefile("data/curated_data/EColi/processed_files/gene_product_molar_mass.json"))
#     reaction_isozymes = Dict(
#         k => Dict(
#             "isozyme_$i" => X.Isozyme(d["gene_product_stoichiometry"], d["kcat_forward"], d["kcat_reverse"])
#             for (i, d) in enumerate(v)
#         )
#         for (k, v) in JSON.parsefile("data/curated_data/EColi/processed_files/isozymes.json")
#     )
#     reaction_isozymes = Dict(
#         rid => Dict("iso" => isozyme[argmax(
#             Dict(
#                 i => iso.kcat_forward /
#                     sum([count * gene_product_molar_masses[gene] for (gene, count) in iso.gene_product_stoichiometry])
#                 for (i, iso) in isozyme
#             )
#         )])
#         for (rid, isozyme) in reaction_isozymes
#     )

#     capacity = JSON.parsefile("data/curated_data/EColi/processed_files/capacity.json")
#     capacity = [(capacity[1][1], string.(capacity[1][2]), capacity[1][3]), (capacity[2][1], string.(capacity[2][2]), capacity[2][3])]
#     c = X.enzyme_constrained_flux_balance_constraints(
#         model;
#         reaction_isozymes,
#         gene_product_molar_masses,
#         capacity,
#     )
# end


# sol = X.optimized_values(
#     c,
#     objective=c.objective.value,
#     optimizer=HiGHS.Optimizer,
# )

# objective_bound = X.relative_tolerance_bound(1)
# objective = c.objective.value

# fva = X.constraints_variability(
#     c * :objective_bound^C.Constraint(objective, objective_bound(sol.objective)),
#     c.fluxes;
#     optimizer=HiGHS.Optimizer,
#     workers = workers(),
# )

# Dict(r => (lb,sol.fluxes[r]) for (r,(lb,ub)) in fva if abs(lb-sol.fluxes[r])>1e-10)

##########################

# Run fva on yeast models 
@everywhere begin
    import AbstractFBCModels as A
    import AbstractFBCModels.CanonicalModel as CM
    import JSONFBCModels as JFBC
    import COBREXA as X
    using JSON, HiGHS
    import ConstraintTrees as C
    import DifferentiableMetabolism as D
    using ElementaryFluxModes
    using ConstraintTrees
    using Random

    idxs = rand(1:343,10)

    organisms = [String(split(x, ".json")[1]) for x in readdir("data/data/fungi343species/models") if !occursin("panmodel",x)][idxs]

    !isdir("data/data/fungi343species/fva") && mkdir("data/data/fungi343species/fva")

    ### calculate sensitivities of all models with a fixed capacity bound of 35% protein/gDW cell
    mass_bound = 350.0
    optimizer = HiGHS.Optimizer
    gene_zero_tol = 1e-8
    flux_zero_tol = 1e-8

    constrs = Dict{String,ConstraintTrees.Tree{ConstraintTrees.Constraint}}()

    for organism in organisms
        # if already done, skip
        "data/data/fungi343species/fva/$organism.json" ∈ readdir("data/data/fungi343species/fva/") && continue

        model = convert(CM.Model, X.load_model("data/data/fungi343species/models/$organism.json"))
        reaction_isozymes = Dict(
            k => Dict(
                "$(k)_$i" => X.Isozyme(
                    d["gene_product_stoichiometry"],
                    d["kcat_forward"],
                    d["kcat_reverse"],
                ) for (i, d) in enumerate(v)
            ) for (k, v) in JSON.parsefile("data/data/fungi343species/isozymes/$organism.json")
        )
        gene_product_molar_masses =
            Dict(x => y for (x, y) in JSON.parsefile("data/data/fungi343species/molar_mass/$organism.json"))
        capacity = [("uncategorized", A.genes(model), mass_bound)]
        # open exchange reactions
        for (r, rxn) in model.reactions
            if contains(rxn.name, "exchange")
                model.reactions[r].lower_bound = rxn.lower_bound * 1000.0
                model.reactions[r].upper_bound = rxn.upper_bound * 1000.0
            end
        end
    
        constrs[organism] = X.enzyme_constrained_flux_balance_constraints(
            model;
            reaction_isozymes,
            gene_product_molar_masses,
            capacity,
        )
    end    

end

for (organism,c) in constrs
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
    isnothing(fva) && continue
    open("data/data/fungi343species/fva/$organism.json","w") do io 
        JSON.print(io,fva)
    end
end


f = JSON.parsefile("data/data/fungi343species/fva/Candida_glabrata.json")

Dict(x=>y for (x,y) in f if all(z->!isnothing(z),y) && abs(y[2]-y[1])>1e-3)
