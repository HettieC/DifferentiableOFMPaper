""" 
Find reactions involving biomass metabolites
"""
function find_reactions_using_metabolite(model,metabolite_id)
    rxns = String[]
    for (r,rxn) in model.reactions
        if haskey(rxn.stoichiometry,metabolite_id)
            push!(rxns,r)
        end
    end
    return rxns
end

"""
Get the bigg database links of a metabolite m.
"""
function get_met_db_bigg(m::Union{String,String31})
    req = nothing
    # make sure the website is still valid 
    try
        req = HTTP.request("GET","http://bigg.ucsd.edu/api/v2/universal/metabolites/$m")
    catch
        req = nothing
    end
    if isnothing(req)
        return nothing
    else
        return JSON.parse(String(req.body))["database_links"]
    end
end
