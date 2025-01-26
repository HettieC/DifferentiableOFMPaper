using COBREXA, JSON, HiGHS, DataFrames, CSV
import AbstractFBCModels.CanonicalModel as CM
import AbstractFBCModels as A
import JSONFBCModels, MATFBCModels
import COBREXA as X
import DifferentiableMetabolism as D
import FastDifferentiation as F
const Ex = F.Node

### calculate sensitivities of all models with a fixed capacity bound of 35% protein/gDW cell
mass_bound = 350.0
optimizer = HiGHS.Optimizer
gene_zero_tol = 1e-8
flux_zero_tol = 1e-8

organisms = [String(split(x, ".json")[1]) for x in readdir("data/fungi343species/models") if !occursin("panmodel",x)]

for organism in organisms
    model = convert(CM.Model, load_model("data/fungi343species/models/$organism.json"))
    reaction_isozymes = Dict(
        k => Dict(
            "$(k)_$i" => X.Isozyme(
                d["gene_product_stoichiometry"],
                d["kcat_forward"],
                d["kcat_reverse"],
            ) for (i, d) in enumerate(v)
        ) for (k, v) in JSON.parsefile("data/fungi343species/isozymes/$organism.json")
    )
    gene_product_molar_masses =
        Dict(x => y for (x, y) in JSON.parsefile("data/fungi343species/molar_mass/$organism.json"))
    capacity = [("uncategorized", A.genes(model), mass_bound)]
    # open exchange reactions
    for (r, rxn) in model.reactions
        if contains(rxn.name, "exchange")
            model.reactions[r].lower_bound = rxn.lower_bound * 1000.0
            model.reactions[r].upper_bound = rxn.upper_bound * 1000.0
        end
    end

    ec_solution = X.enzyme_constrained_flux_balance_analysis(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
        optimizer=HiGHS.Optimizer
    )
    

    #if capacity bound not hit, skip this model
    if abs(ec_solution.gene_product_capacity.uncategorized - mass_bound) > 1e-5
        println(organism)
        break
    end

    # prune the optimal solution
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

    parameter_isozymes = Dict{String,Dict{String,X.IsozymeT{Ex}}}()
    parameter_values = Dict{Symbol,Float64}()
    for (rid,iso) in pruned_reaction_isozymes
        parameter_values[Symbol(rid)] = collect(values(iso))[1].kcat_forward
        parameter_isozymes[rid] = Dict(
            "isozyme" => X.IsozymeT{Ex}(
                b.gene_product_stoichiometry, 
                Ex(Symbol(rid)),
                nothing
            )
            for (_,b) in iso
        )
    end
    
    rid_kcat = Dict(k => Ex(Symbol(k)) for (k,_) in parameter_values)

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
    );
    df = DataFrame(sens, :auto);
    rename!(df, string.(parameters));
    insertcols!(df, 1, :vid => vids);
    CSV.write("data/fungi343species/fixed_EC_sensitivities/$organism.csv", df);
end
