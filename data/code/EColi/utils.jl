function get_protein_stoichiometry!(model, complex, proteome_data)

    #: protein stoich map, infer from uniprot
    mer_map = [
        "Homotetramer" 4
        "Homodimer" 2
        "Homotrimer" 3
        "Homohexamer" 6
        "Homopentamer" 5
        "Homodecamer" 10
        "Homooctamer" 8
        "Homoheptamer" 7
        "Homododecamer" 12
        "Homomonomer" 1
        "Monomer" 1
    ]

    #: infer protein stoichiometry from uniprot annotations
    protein_stoichiometry = Dict()
    for rid in reactions(model)
        if !isnothing(A.reaction_gene_association_dnf(model, rid))
            grrs = A.reaction_gene_association_dnf(model, rid)
            smer = []
            for grr in grrs
                if length(grr) == 1 # only assign homomers
                    gid = grr[1]
                    if haskey(proteome_data, gid)
                        mer = proteome_data[gid][2]
                        if any(startswith.(Ref(mer), mer_map[:, 1]))
                            idx = findfirst(x -> startswith(mer, x), mer_map[:, 1])
                            push!(smer, [mer_map[idx, 2]])
                        else
                            push!(smer, [1.0])
                        end
                    else # no data
                        push!(smer, [1.0])
                    end
                else # assume complexes have uni-stoichiometry, manually fix later
                    push!(smer, fill(1.0, length(grr)))
                end
            end
            isempty(smer[1]) && continue
            protein_stoichiometry[rid] = smer
        end
    end

    #:fix complex stoichiometry, use ComplexPortal database

    for rid in reactions(model)
        isnothing(A.reaction_gene_association_dnf(model, rid)) && continue
        grrs = A.reaction_gene_association_dnf(model, rid)
        length(first(grrs)) == 1 && continue # skip monomers

        accurate_complex_idxs = Int[]
        for (i, grr) in enumerate(grrs)
            stoichs = []
            for v in values(complex)
                if length(intersect(collect(keys(v)), grr)) == length(grr) &&
                   length(intersect(collect(keys(v)), grr)) == length(v)
                    push!(stoichs, v)
                end
            end
            if length(stoichs) > 1
                @warn("Uncertain which complex to choose for reaction", rid)
            elseif length(stoichs) == 1
                stoich = first(stoichs)
                protein_stoichiometry[rid][i] = [
                    get(stoich, gid, 1.0) == 0 ? 1.0 : get(stoich, gid, 1.0) for gid in grr
                ]
                push!(accurate_complex_idxs, i)
            end
        end

        #: remove grrs for complexes not found in the database (conservative)
        if !isempty(accurate_complex_idxs)
            rem_idxs = filter(x -> x ∉ accurate_complex_idxs, 1:length(grrs))
            deleteat!(model.reactions[rid].gene_association_dnf, rem_idxs)
            deleteat!(protein_stoichiometry[rid], rem_idxs)
        end
    end

    return protein_stoichiometry
end
