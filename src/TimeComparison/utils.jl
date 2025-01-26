"""
Calculate EFMs and the change in lambda
"""
function calculate_new_lambda(pruned_model, N, pruned_solution, atpm, biomass)
    atpm_idx = findfirst(x -> x == atpm, A.reactions(pruned_model))
    biomass_idx = findfirst(x -> x == biomass, A.reactions(pruned_model))
    fixed_fluxes = [atpm_idx, biomass_idx]
    flux_values = [pruned_solution.fluxes[atpm], pruned_solution.fluxes[biomass]]

    OFMs = get_ofms(Matrix(N), fixed_fluxes, flux_values)
    OFM_dicts = [
        Dict(A.reactions(pruned_model) .=> OFMs[:, 1]),
        Dict(A.reactions(pruned_model) .=> OFMs[:, 2])
    ]
    # scale to have one unit flux through biomass
    OFM_dicts = [
        Dict(x => y / OFM_dicts[1][biomass] for (x, y) in OFM_dicts[1]),
        Dict(x => y / OFM_dicts[2][biomass] for (x, y) in OFM_dicts[2])
    ]
    ## calculate lambda
    rxn1 = [x for (x, y) in OFM_dicts[1] if OFM_dicts[2][x] == 0 && y > 1e-5][1]
    rxn2 = [x for (x, y) in OFM_dicts[2] if OFM_dicts[1][x] == 0 && y > 1e-5][1]
    M = [
        OFM_dicts[1][rxn1] OFM_dicts[2][rxn1];
        OFM_dicts[1][rxn2] OFM_dicts[2][rxn2]
    ]

    v = [
        pruned_solution.fluxes[rxn1];
        pruned_solution.fluxes[rxn2]
    ]

    return M \ v
end

"""
Get the timing of doing finite diff to get sensitivity of efms to all kcats 
"""
function finite_diff_ofms(mod, isos, mass, cap, atpm, biomass, glc_exch)
    model = convert(CM.Model, X.load_model(mod))
    model.reactions[glc_exch].lower_bound = -1000.0
    gene_product_molar_masses = Dict(x => y for (x, y) in JSON.parsefile(mass))
    reaction_isozymes = Dict(
        k => Dict(
            "isozyme_$i" => X.Isozyme(d["gene_product_stoichiometry"], d["kcat_forward"], d["kcat_reverse"])
            for (i, d) in enumerate(v)
        )
        for (k, v) in JSON.parsefile(isos)
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
    capacity = JSON.parsefile(cap)
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
    )

    N = A.stoichiometry(pruned_model)

    ## loop over parameters and change kcats one by one
    delta_pos = 1.001
    delta_neg = 0.999

    d_lambda = Dict{String,Vector{Float64}}()
    d_vr = Dict{String,Float64}()

    timing = @benchmark begin
        i = 0
        for r in collect(keys($pruned_reaction_isozymes))
            i += 1
            $pruned_reaction_isozymes[r]["iso"].kcat_forward *= $delta_pos
            # solve the gecko model with new isozyme
            ec_solution_new = X.enzyme_constrained_flux_balance_analysis(
                $pruned_model;
                reaction_isozymes=$pruned_reaction_isozymes,
                gene_product_molar_masses=$gene_product_molar_masses,
                capacity=$capacity,
                optimizer=$HiGHS.Optimizer
            )

            vr_pos = ec_solution_new.fluxes[$biomass]

            # calculate new EFMs based on the new pgm 
            flux_values = [ec_solution_new.fluxes[$atpm], ec_solution_new.fluxes[$biomass]]

            lambda_new_pos = calculate_new_lambda($pruned_model, $N, ec_solution_new, $atpm, $biomass)

            $pruned_reaction_isozymes[r]["iso"].kcat_forward *= $delta_neg / $delta_pos

            # solve the gecko model with new isozyme
            ec_solution_new = X.enzyme_constrained_flux_balance_analysis(
                $pruned_model;
                reaction_isozymes=$pruned_reaction_isozymes,
                gene_product_molar_masses=$gene_product_molar_masses,
                capacity=$capacity,
                optimizer=HiGHS.Optimizer
            )

            vr_neg = ec_solution_new.fluxes[$biomass]

            # calculate new EFMs based on the new pgm 
            flux_values = [ec_solution_new.fluxes[$atpm], ec_solution_new.fluxes[$biomass]]

            lambda_new_neg = calculate_new_lambda($pruned_model, $N, ec_solution_new, $atpm, $biomass)


            #return to original value
            $pruned_reaction_isozymes[r]["iso"].kcat_forward /= $delta_neg

            $d_lambda[r] = (lambda_new_pos .- lambda_new_neg) ./ (($delta_pos - $delta_neg) * $pruned_reaction_isozymes[r]["iso"].kcat_forward)
            $d_vr[r] = (
                vr_pos - vr_neg
            ) / (($delta_pos - $delta_neg) * $pruned_reaction_isozymes[r]["iso"].kcat_forward)
        end
    end
    return timing
end

"""
Get the timing of DiffMet to get efm sensitivity
"""
function diff_met_ofm(mod, isos, mass, cap, atpm, biomass, glc_exch)
    model = convert(CM.Model, X.load_model(mod))
    model.reactions[glc_exch].lower_bound = -1000.0
    gene_product_molar_masses = Dict(x => y for (x, y) in JSON.parsefile(mass))
    reaction_isozymes = Dict(
        k => Dict(
            "isozyme_$i" => X.Isozyme(d["gene_product_stoichiometry"], d["kcat_forward"], d["kcat_reverse"])
            for (i, d) in enumerate(v)
        )
        for (k, v) in JSON.parsefile(isos)
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
    capacity = JSON.parsefile(cap)
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
    )

    ec_solution_pruned = X.enzyme_constrained_flux_balance_analysis(
        pruned_model;
        reaction_isozymes=pruned_reaction_isozymes,
        gene_product_molar_masses,
        capacity,
        optimizer=HiGHS.Optimizer
    )

    # calculate OFMs
    parameter_values = Dict(Symbol(x) => y["iso"].kcat_forward for (x, y) in pruned_reaction_isozymes)
    rid_kcat = Dict(k => Ex(Symbol(k)) for (k, _) in parameter_values)
    parameter_isozymes = Dict(
        x => Dict(
            "iso" => X.IsozymeT{Ex}(
                y["iso"].gene_product_stoichiometry,
                rid_kcat[Symbol(x)],
                nothing
            )
        )
        for (x, y) in pruned_reaction_isozymes
    )

    pkm = X.enzyme_constrained_flux_balance_constraints( # kinetic model
        pruned_model;
        reaction_isozymes=parameter_isozymes,
        gene_product_molar_masses,
        capacity,
    )
    pruned_solution = D.optimized_values(
        pkm,
        parameter_values;
        objective=pkm.objective.value,
        optimizer=HiGHS.Optimizer,
    )

    # calculate EFMs
    N = A.stoichiometry(pruned_model)

    # atpm, biomass
    atpm_idx = findfirst(x -> x == atpm, A.reactions(pruned_model))
    biomass_idx = findfirst(x -> x == biomass, A.reactions(pruned_model))
    fixed_fluxes = [atpm_idx, biomass_idx]
    flux_values = [pruned_solution.tree.fluxes[atpm], pruned_solution.tree.fluxes[biomass]]

    OFMs = get_ofms(Matrix(N), fixed_fluxes, flux_values)

    OFM_dicts = [
        Dict(A.reactions(pruned_model) .=> OFMs[:, 1]),
        Dict(A.reactions(pruned_model) .=> OFMs[:, 2])
    ]
    # scale to have one unit flux through biomass
    OFM_dicts = [
        Dict(x => y / OFM_dicts[1][biomass] for (x, y) in OFM_dicts[1]),
        Dict(x => y / OFM_dicts[2][biomass] for (x, y) in OFM_dicts[2])
    ]

    parameter_values = Dict(Symbol(x) => first(values(y)).kcat_forward for (x, y) in pruned_reaction_isozymes)
    parameters = Ex.(collect(keys(parameter_values)))
    rid_pid = Dict(rid => [iso.kcat_forward for (k, iso) in v][1] for (rid, v) in parameter_isozymes)
    rid_gcounts = Dict(rid => [v.gene_product_stoichiometry for (k, v) in d][1] for (rid, d) in pruned_reaction_isozymes)

    sens_efm = differentiate_efm(OFM_dicts, parameters, rid_pid, parameter_values, rid_gcounts, capacity, gene_product_molar_masses, HiGHS.Optimizer)

    timing = @benchmark differentiate_efm($OFM_dicts, $parameters, $rid_pid, $parameter_values, $rid_gcounts, $capacity, $gene_product_molar_masses, $HiGHS.Optimizer) samples = 10

    return timing
end

"""
DiffMet timing of model. Return: whole diffmet time, GD time, number of reactions in solution
"""
function diff_met(mod::String, isos::String, mass::String, cap::String,glc_exch)
    model = convert(CM.Model, X.load_model(mod))
    model.reactions[glc_exch].lower_bound = -1000.0
    gene_product_molar_masses = Dict(x => y for (x, y) in JSON.parsefile(mass))
    reaction_isozymes = Dict(
        k => Dict(
            "isozyme_$i" => X.Isozyme(d["gene_product_stoichiometry"], d["kcat_forward"], d["kcat_reverse"])
            for (i, d) in enumerate(v)
        )
        for (k, v) in JSON.parsefile(isos)
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
    capacity = JSON.parsefile(cap)
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
    )
    parameter_values = Dict(Symbol(x) => y["iso"].kcat_forward for (x, y) in pruned_reaction_isozymes)
    rid_kcat = Dict(k => Ex(Symbol(k)) for (k, _) in parameter_values)
    parameter_isozymes = Dict(
        x => Dict(
            "iso" => X.IsozymeT{Ex}(
                y["iso"].gene_product_stoichiometry,
                rid_kcat[Symbol(x)],
                nothing
            )
        )
        for (x, y) in pruned_reaction_isozymes
    )

    pkm = X.enzyme_constrained_flux_balance_constraints( # kinetic model
        pruned_model;
        reaction_isozymes=parameter_isozymes,
        gene_product_molar_masses,
        capacity,
    )
    pruned_solution = D.optimized_values(
        pkm,
        parameter_values;
        objective=pkm.objective.value,
        optimizer=HiGHS.Optimizer,
    )

    parameters = collect(keys(parameter_values))


    # do whole diff met
    dm = @benchmark begin
        pkm_kkt, vids = D.differentiate_prepare_kkt($pkm, $pkm.objective.value, $parameters)

        sens = D.differentiate_solution(
            pkm_kkt,
            $pruned_solution.primal_values,
            $pruned_solution.equality_dual_values,
            $pruned_solution.inequality_dual_values,
            $parameter_values,
            scale=true, # unitless sensitivities
        )
    end
    return dm, length(pruned_reaction_isozymes)
end

"""
Use finite difference on a model.
"""
function finite_difference(mod, isos, mass, cap, glc_exch)
    model = convert(CM.Model, X.load_model(mod))
    model.reactions[glc_exch].lower_bound = -1000.0
    gene_product_molar_masses = Dict(x => y for (x, y) in JSON.parsefile(mass))
    reaction_isozymes = Dict(
        k => Dict(
            "isozyme_$i" => X.Isozyme(d["gene_product_stoichiometry"], d["kcat_forward"], d["kcat_reverse"])
            for (i, d) in enumerate(v)
        )
        for (k, v) in JSON.parsefile(isos)
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
    capacity = JSON.parsefile(cap)
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
    )


    ## loop over parameters and change kcats one by one
    delta_pos = 1.001
    delta_neg = 0.999

    d_vr = Dict{String,Float64}()

    timing = @benchmark begin
        i = 0
        for r in collect(keys($pruned_reaction_isozymes))
            i += 1
            $pruned_reaction_isozymes[r]["iso"].kcat_forward *= $delta_pos
            # solve the gecko model with new isozyme
            ec_solution_new = X.enzyme_constrained_flux_balance_analysis(
                $pruned_model;
                reaction_isozymes=$pruned_reaction_isozymes,
                gene_product_molar_masses=$gene_product_molar_masses,
                capacity=$capacity,
                optimizer=$HiGHS.Optimizer
            )

            vr_pos = ec_solution_new.objective

            $pruned_reaction_isozymes[r]["iso"].kcat_forward *= $delta_neg / $delta_pos

            # solve the gecko model with new isozyme
            ec_solution_new = X.enzyme_constrained_flux_balance_analysis(
                $pruned_model;
                reaction_isozymes=$pruned_reaction_isozymes,
                gene_product_molar_masses=$gene_product_molar_masses,
                capacity=$capacity,
                optimizer=HiGHS.Optimizer
            )

            vr_neg = ec_solution_new.objective

            #return to original value
            $pruned_reaction_isozymes[r]["iso"].kcat_forward /= $delta_neg

            $d_vr[r] = (
                vr_pos - vr_neg
            ) / (($delta_pos - $delta_neg) * $pruned_reaction_isozymes[r]["iso"].kcat_forward)
        end
    end
    return timing
end
