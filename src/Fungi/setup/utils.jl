using FastaIO
function get_gene_product_molar_mass(organism,model)
    # molecular weight of the amino acids in g/mol (1 Da = 1 g/mol)
    AAmw = Dict(
        "A" => 89.1,
        "R" => 174.2,
        "N" => 132.1,
        "D" => 133.1,
        "C" => 121.2,
        "E" => 147.1,
        "Q" => 146.2,
        "G" => 75.1,
        "H" => 155.2,
        "I" => 131.2,
        "L" => 131.2,
        "K" => 146.2,
        "M" => 149.2,
        "F" => 165.2,
        "P" => 115.1,
        "S" => 105.1,
        "T" => 119.1,
        "W" => 204.2,
        "Y" => 181.2,
        "V" => 117.1,
        "X" => 136.9,
    )

    # load the protein sequences from fasta file
    fasta = readfasta("data/fungi343species/Proteinfasta/$organism.fasta")
    gene_product_molar_mass = Dict{String,Float64}()
    for (name, seq) in fasta
        mw = 0
        for AA in seq
            mw += AAmw["$AA"]
        end
        gene_product_molar_mass[name] = mw / 1000.0
    end
    avg_mw = sum(values(gene_product_molar_mass)) / length(gene_product_molar_mass)
    for (g,gene) in model.genes
        if !haskey(gene_product_molar_mass,g)
            gene_product_molar_mass[g] = avg_mw
        end
    end
    return gene_product_molar_mass
end

"""
Find a mass bound that gets hit.
"""
function find_mass_bound(model, mass_bound, growth, reaction_parameter_isozymes, rxn_isozymes, gene_product_molar_masses, capacity, optimizer)
    parameter_values = Dict{Symbol,Float64}()
    for (r, iso) in rxn_isozymes
        for (k, v) in iso
            parameter_values[Symbol(k, :_f)] = v.kcat_forward
            parameter_values[Symbol(k, :_r)] = v.kcat_reverse
        end
    end
    variable_values = Dict(Symbolics.variable(k) => v for (k, v) in parameter_values);
    km = build_kinetic_model(
        model;
        reaction_isozymes = reaction_parameter_isozymes,
        gene_product_molar_masses,
        capacity,
    )
    ec_solution, _, _, _ = optimized_constraints_with_parameters(
        km,
        variable_values;
        objective = km.objective.value,
        optimizer,
    )

    if 0.95 < ec_solution.objective/growth < 1.05
        return mass_bound, ec_solution
    else
        new_mass_bound = mass_bound*growth/ec_solution.objective
        capacity = [(capacity[1][1],capacity[1][2],new_mass_bound)]
        find_mass_bound(model, new_mass_bound, growth, reaction_parameter_isozymes, rxn_isozymes, gene_product_molar_masses, capacity, optimizer)
    end
end

function finite_diff(pruned_model, pruned_reaction_isozymes, pruned_ec_solution, gene_product_molar_masses, capacity, optimizer)
    sens = Matrix(undef,length(pruned_model.reactions),length(pruned_reaction_isozymes))
    pids = String[]
    i = 0 
    for (k,v) in pruned_reaction_isozymes
        for (x,y) in v
            i +=1
            pruned_reaction_isozymes[k][x] = Isozyme(
                y.gene_product_stoichiometry,
                ec_solution.fluxes[k] > 0 ? 1.01*y.kcat_forward : 0.0,
                ec_solution.fluxes[k] < 0 ? 1.01*y.kcat_reverse : 0.0
            )

            ec_solution_new = enzyme_constrained_flux_balance_analysis(
                pruned_model;
                reaction_isozymes = pruned_reaction_isozymes,
                gene_product_molar_masses,
                capacity,
                optimizer,
            )
            pruned_reaction_isozymes[k][x] = Isozyme(
                y.gene_product_stoichiometry,
                y.kcat_forward,
                y.kcat_reverse
            )
            push!(pids,x)
            sens[:,i] = (collect(values(pruned_ec_solution.fluxes)) .- collect(values(ec_solution_new.fluxes)))/ (ec_solution.fluxes[k] > 0 ? 0.01*y.kcat_forward : 0.01*y.kcat_reverse)
        end
    end
    return sens, pids
end

"""
TO DO Check if reactions have atom mass balance
"""
function atom_balances(model)
    balances = Dict{String,Bool}()
    for (r,rxn) in model.reactions
        if length(rxn.stoichiometry) == 1 
            balances[r] = false
        else
            balance = Dict{String,Union{Int,Float64}}()
            for (m,s) in rxn.stoichiometry
                for (x,y) in model.metabolites[m].formula
                    if haskey(balance,x)
                        balance[x] += y*s
                    else
                        balance[x] = y*s
                    end
                end
            end
            if all(((k,v),) -> v==0, balance)
                return true, balance
            else  
                return false, balance
            end
        end
    end
    return balances
end
