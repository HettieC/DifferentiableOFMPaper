using FastaIO

function get_gene_product_molar_mass(model, path_to_fasta)
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
    fasta = readfasta(path_to_fasta)
    gene_product_molar_mass = Dict{String,Float64}()
    for (name, seq) in fasta
        mw = 0
        for AA in seq
            mw += AAmw["$AA"]
        end
        gene_product_molar_mass[name] = mw / 1000.0
    end
    avg_mw = sum(values(gene_product_molar_mass)) / length(gene_product_molar_mass)
    for (g, gene) in model.genes
        if !haskey(gene_product_molar_mass, g)
            gene_product_molar_mass[g] = avg_mw
        end
    end
    return gene_product_molar_mass
end
