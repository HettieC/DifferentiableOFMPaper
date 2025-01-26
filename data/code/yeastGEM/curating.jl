#= 
Using the Chalmers Yeast9 model
Zhang, C. et al. Yeast9: a consensus yeast metabolic model enables quantitative analysis of 
cellular metabolism by incorporating big data. bioRxiv (2023) doi:10.1101/2023.12.03.569754
version 9.0.0 
=#

using COBREXA, MATFBCModels
import AbstractFBCModels.CanonicalModel as CM
import AbstractFBCModels as A
using DataFrames, CSV, JSONFBCModels
include("utils.jl")
#### download like in cobrexa docs
download_model(
    "https://github.com/SysBioChalmers/yeast-GEM/archive/refs/tags/v9.0.1.zip",
    "data/data/yeastGEM/yeast-GEM.mat",
    "7c4643468f18442d78eea81d68343eb7ff814b5eac91c3cbeb67c7b736c5f640"
)

model = convert(CM.Model,load_model("data/data/yeastGEM/yeast-GEM.mat"))

biomass_rxns = String[]
for (m,s) in model.reactions["r_4041"].stoichiometry 
    if s < 0 && m ∉ ["s_0803","s_0434"]
        append!(biomass_rxns,find_reactions_using_metabolite(model,m))
    end
end
for r in biomass_rxns
    for (k,v) in model.reactions[r].stoichiometry 
        if abs(v) < 1e-4
            model.reactions[r].stoichiometry[k] = 0.0
        end
    end
end

# add smiles annotations
smiles_df = DataFrame(CSV.File("data/data/yeastGEM/databases/smilesDB.tsv",header=false))
smiles_dict = Dict(Pair.(smiles_df.Column1,smiles_df.Column2))
mets_smiles = Dict(x => smiles_dict[y.name] for (x,y) in model.metabolites)
for (k,v) in mets_smiles
    if !ismissing(v)
        model.metabolites[k].annotations["smiles"] = [v]
    end
end


### get smiles/inchi/kegg for metabolites without smiles in the smilesDB 
bigg_dict_first = CSV.File("data/data/yeastGEM/databases/BiGGmetDictionary.csv",header=false) |> Dict
bigg_dict_second = CSV.File("data/data/yeastGEM/databases/BiGGmetDictionary_newIDs.csv",header=false) |> Dict
bigg_dict = merge(
    Dict("M_$(split(x,"[")[1])" => y for (x,y) in bigg_dict_first),
    Dict("M_$(split(x,"[")[1])" => y for (x,y) in bigg_dict_second)
)
mets_to_get = Dict(x => bigg_dict[x] for (x,y) in mets_smiles if haskey(bigg_dict,x) && ismissing(y))
mets_no_bigg = [x for (x,y) in mets_smiles if !haskey(bigg_dict,x) && ismissing(y)]
for (m,bigg) in mets_to_get
    db = get_met_db_bigg(bigg)
    isnothing(db) && continue 
    if haskey(db,"KEGG Compound")
        mets_smiles[m] = db["KEGG Compound"][1]["id"]
        model.metabolites[m].annotations["KEGG Compound"] = [db["KEGG Compound"][1]["id"],db["KEGG Compound"][1]["link"]]
    elseif haskey(db,"MetaNetX (MNX) Chemical")
        mets_smiles[m] = db["MetaNetX (MNX) Chemical"][1]["link"]
        model.metabolites[m].annotations["MetaNetX (MNX) Chemical"] = [db["MetaNetX (MNX) Chemical"][1]["id"],db["MetaNetX (MNX) Chemical"][1]["link"]]
    end
end

### manually add the inchi strings where needed
mets_smiles["s_3424"] = "InChI=1S/C21H44NO7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-21(24)27-18-20(23)19-29-30(25,26)28-17-16-22/h20,23H,2-19,22H2,1H3,(H,25,26)/t20-/m1/s1"
mets_smiles["s_3520"] = "InChI=1S/C39H72O5/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-31-33-38(41)43-36-37(35-40)44-39(42)34-32-30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,37,40H,3-12,17-36H2,1-2H3"
mets_smiles["s_2982"] = "InChI=1S/C35H64O5/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-34(37)39-32-33(31-36)40-35(38)30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,33,36H,3-12,17-32H2,1-2H3"
mets_smiles["s_3632"] = "InChI=1S/C19H39O7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-19(21)25-16-18(20)17-26-27(22,23)24/h18,20H,2-17H2,1H3,(H2,22,23,24)/p-2/t18-/m1/s1"
mets_smiles["s_3687"] = "InChI=1S/C21H42NO7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-21(24)27-18-20(23)19-29-30(25,26)28-17-16-22/h9-10,20,23H,2-8,11-19,22H2,1H3,(H,25,26)/t20-/m1/s1"
mets_smiles["s_3648"] = "InChI=1S/C19H39O7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-19(21)25-16-18(20)17-26-27(22,23)24/h18,20H,2-17H2,1H3,(H2,22,23,24)/p-2/t18-/m1/s1"
mets_smiles["s_3706"] = "InChI=1S/C37H70NO8P/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-36(39)43-33-35(34-45-47(41,42)44-32-31-38)46-37(40)30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,35H,3-12,17-34,38H2,1-2H3,(H,41,42)"
mets_smiles["s_3010"] = "InChI=1S/C39H72O5/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-31-33-38(41)43-36-37(35-40)44-39(42)34-32-30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,37,40H,3-12,17-36H2,1-2H3"
mets_smiles["s_3241"] = "InChI=1S/C73H134O17P2/c1-5-9-13-17-21-25-29-33-37-41-45-49-53-57-70(75)83-63-68(89-72(77)59-55-51-47-43-39-35-31-27-23-19-15-11-7-3)65-87-91(79,80)85-61-67(74)62-86-92(81,82)88-66-69(90-73(78)60-56-52-48-44-40-36-32-28-24-20-16-12-8-4)64-84-71(76)58-54-50-46-42-38-34-30-26-22-18-14-10-6-2/h25-32,67-69,74H,5-24,33-66H2,1-4H3,(H,79,80)(H,81,82)/p-2"
mets_smiles["s_3150"] = "InChI=1S/C37H70NO8P/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-36(39)43-33-35(34-45-47(41,42)44-32-31-38)46-37(40)30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,35H,3-12,17-34,38H2,1-2H3,(H,41,42)"
mets_smiles["s_3047"] = "InChI=1S/C35H64O5/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-34(37)39-32-33(31-36)40-35(38)30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,33,36H,3-12,17-32H2,1-2H3"
mets_smiles["s_2935"] = "InChI=1S/C19H39O7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-19(21)25-16-18(20)17-26-27(22,23)24/h18,20H,2-17H2,1H3,(H2,22,23,24)/p-2/t18-/m1/s1"
mets_smiles["s_3636"] = "InChI=1S/C21H43O7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-21(23)27-18-20(22)19-28-29(24,25)26/h20,22H,2-19H2,1H3,(H2,24,25,26)/p-2/t20-/m1/s1"
mets_smiles["s_3500"] = "InChI=1S/C35H64O5/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-34(37)39-32-33(31-36)40-35(38)30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,33,36H,3-12,17-32H2,1-2H3"
mets_smiles["s_0302"] = "InChI=1S/C8H16N3O8P/c9-5(1-10-3-12)11-8-7(14)6(13)4(19-8)2-18-20(15,16)17/h3-4,6-8,13-14H,1-2H2,(H2,9,11)(H,10,12)(H2,15,16,17)/p-1/t4-,6-,7-,8-/m1/s1"
mets_smiles["s_3512"] = "InChI=1S/C39H72O5/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-31-33-38(41)43-36-37(35-40)44-39(42)34-32-30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,37,40H,3-12,17-36H2,1-2H3"
mets_smiles["s_3640"] = "InChI=1S/C19H39O7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-19(21)25-16-18(20)17-26-27(22,23)24/h18,20H,2-17H2,1H3,(H2,22,23,24)/p-2/t18-/m1/s1"
mets_smiles["s_3185"] = "InChI=1S/C37H70NO8P/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-36(39)43-33-35(34-45-47(41,42)44-32-31-38)46-37(40)30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,35H,3-12,17-34,38H2,1-2H3,(H,41,42)"
mets_smiles["s_2974"] = "InChI=1S/C39H72O5/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-31-33-38(41)43-36-37(35-40)44-39(42)34-32-30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,37,40H,3-12,17-36H2,1-2H3"
mets_smiles["s_3057"] = "InChI=1S/C39H72O5/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-31-33-38(41)43-36-37(35-40)44-39(42)34-32-30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,37,40H,3-12,17-36H2,1-2H3"
mets_smiles["s_2891"] = "InChI=1S/C45H80N7O17P3S/c1-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18-19-20-21-22-23-24-25-36(54)73-29-28-47-35(53)26-27-48-43(57)40(56)45(2,3)31-66-72(63,64)69-71(61,62)65-30-34-39(68-70(58,59)60)38(55)44(67-34)52-33-51-37-41(46)49-32-50-42(37)52/h24-25,32-34,38-40,44,55-56H,4-23,26-31H2,1-3H3,(H,47,53)(H,48,57)(H,61,62)(H,63,64)(H2,46,49,50)(H2,58,59,60)/p-4/b25-24+/t34?,38?,39?,40?,44-/m0/s1"
mets_smiles["s_3168"] = "InChI=1S/C37H70NO8P/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-36(39)43-33-35(34-45-47(41,42)44-32-31-38)46-37(40)30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,35H,3-12,17-34,38H2,1-2H3,(H,41,42)"
mets_smiles["s_2956"] = "	InChI=1S/C35H65O8P/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-34(36)41-31-33(32-42-44(38,39)40)43-35(37)30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,33H,3-12,17-32H2,1-2H3,(H2,38,39,40)/p-1"
mets_smiles["s_0404"] = "CC([NH3+])C(=O)OC1C(O)C(C[*])OC1n1cnc2c(N)ncnc21"
mets_smiles["s_0670"] = "InChI=1S/C34H54O6/c1-19(2)20(3)7-8-21(4)25-11-12-26-24-10-9-22-17-23(13-15-33(22,5)27(24)14-16-34(25,26)6)39-32-31(38)30(37)29(36)28(18-35)40-32/h7-10,19-21,23,25-32,35-38H,11-18H2,1-6H3/b8-7+/t20-,21+,23-,25+,26-,27-,28+,29+,30-,31+,32+,33-,34+/m0/s1"
mets_smiles["s_3514"] = "InChI=1S/C35H64O5/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-34(37)39-32-33(31-36)40-35(38)30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,33,36H,3-12,17-32H2,1-2H3"
mets_smiles["s_2944"] = "InChI=1S/C19H39O7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-19(21)25-16-18(20)17-26-27(22,23)24/h18,20H,2-17H2,1H3,(H2,22,23,24)/p-2/t18-/m1/s1"
mets_smiles["s_3522"] = "InChI=1S/C35H64O5/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-34(37)39-32-33(31-36)40-35(38)30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,33,36H,3-12,17-32H2,1-2H3"
mets_smiles["s_3644"] = "InChI=1S/C21H43O7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-21(23)27-18-20(22)19-28-29(24,25)26/h20,22H,2-19H2,1H3,(H2,24,25,26)/p-2/t20-/m1/s1"
mets_smiles["s_3468"] = "InChI=1S/C21H42NO7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-21(24)27-18-20(23)19-29-30(25,26)28-17-16-22/h9-10,20,23H,2-8,11-19,22H2,1H3,(H,25,26)/t20-/m1/s1"
mets_smiles["s_2814"] = "InChI=1S/C45H80N7O17P3S/c1-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18-19-20-21-22-23-24-25-36(54)73-29-28-47-35(53)26-27-48-43(57)40(56)45(2,3)31-66-72(63,64)69-71(61,62)65-30-34-39(68-70(58,59)60)38(55)44(67-34)52-33-51-37-41(46)49-32-50-42(37)52/h24-25,32-34,38-40,44,55-56H,4-23,26-31H2,1-3H3,(H,47,53)(H,48,57)(H,61,62)(H,63,64)(H2,46,49,50)(H2,58,59,60)/p-4/b25-24+/t34?,38?,39?,40?,44-/m0/s1"
mets_smiles["s_3682"] = "InChI=1S/C39H72O5/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-31-33-38(41)43-36-37(35-40)44-39(42)34-32-30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,37,40H,3-12,17-36H2,1-2H3"
mets_smiles["s_2929"] = "InChI=1S/C33H56N7O17P3S/c1-4-5-6-7-8-9-10-11-12-13-24(42)61-17-16-35-23(41)14-15-36-31(45)28(44)33(2,3)19-54-60(51,52)57-59(49,50)53-18-22-27(56-58(46,47)48)26(43)32(55-22)40-21-39-25-29(34)37-20-38-30(25)40/h11-12,20-22,26-28,32,43-44H,4-10,13-19H2,1-3H3,(H,35,41)(H,36,45)(H,49,50)(H,51,52)(H2,34,37,38)(H2,46,47,48)/p-4/b12-11-/t22-,26-,27-,28+,32-/m1/s1"
mets_smiles["s_3691"] = "InChI=1S/C23H46NO7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-23(26)29-20-22(25)21-31-32(27,28)30-19-18-24/h9-10,22,25H,2-8,11-21,24H2,1H3,(H,27,28)/b10-9+"
mets_smiles["s_3528"] = "InChI=1S/C35H65O8P/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-34(36)41-31-33(32-42-44(38,39)40)43-35(37)30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,33H,3-12,17-32H2,1-2H3,(H2,38,39,40)/p-1"
mets_smiles["s_3000"] = "InChI=1S/C35H64O5/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-34(37)39-32-33(31-36)40-35(38)30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,33,36H,3-12,17-32H2,1-2H3"
mets_smiles["s_2946"] = "InChI=1S/C21H43O7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-21(23)27-18-20(22)19-28-29(24,25)26/h20,22H,2-19H2,1H3,(H2,24,25,26)/p-2/t20-/m1/s1"
mets_smiles["s_3466"] = "InChI=1S/C21H44NO7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-21(24)27-18-20(23)19-29-30(25,26)28-17-16-22/h20,23H,2-19,22H2,1H3,(H,25,26)/t20-/m1/s1"
mets_smiles["s_2992"] = "InChI=1S/C39H72O5/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-31-33-38(41)43-36-37(35-40)44-39(42)34-32-30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,37,40H,3-12,17-36H2,1-2H3"
mets_smiles["s_3132"] = "InChI=1S/C37H70NO8P/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-36(39)43-33-35(34-45-47(41,42)44-32-31-38)46-37(40)30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,35H,3-12,17-34,38H2,1-2H3,(H,41,42)"
mets_smiles["s_2937"] = "InChI=1S/C21H43O7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-21(23)27-18-20(22)19-28-29(24,25)26/h20,22H,2-19H2,1H3,(H2,24,25,26)/p-2/t20-/m1/s1"
mets_smiles["s_3652"] = "InChI=1S/C21H43O7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-21(23)27-18-20(22)19-28-29(24,25)26/h20,22H,2-19H2,1H3,(H2,24,25,26)/p-2/t20-/m1/s1"
mets_smiles["s_2969"] = "InChI=1S/C35H64O5/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-34(37)39-32-33(31-36)40-35(38)30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,33,36H,3-12,17-32H2,1-2H3"
mets_smiles["s_2999"] = "InChI=1S/C35H65O8P/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-34(36)41-31-33(32-42-44(38,39)40)43-35(37)30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,33H,3-12,17-32H2,1-2H3,(H2,38,39,40)/p-1"
mets_smiles["s_2981"] = "InChI=1S/C35H65O8P/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-34(36)41-31-33(32-42-44(38,39)40)43-35(37)30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,33H,3-12,17-32H2,1-2H3,(H2,38,39,40)/p-1"
mets_smiles["s_3685"] = "InChI=1S/C21H44NO7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-21(24)27-18-20(23)19-29-30(25,26)28-17-16-22/h20,23H,2-19,22H2,1H3,(H,25,26)/t20-/m1/s1"
mets_smiles["s_3425"] = "InChI=1S/C21H42NO7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-21(24)27-18-20(23)19-29-30(25,26)28-17-16-22/h9-10,20,23H,2-8,11-19,22H2,1H3,(H,25,26)/t20-/m1/s1"
mets_smiles["s_3467"] = "InChI=1S/C37H70NO8P/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-36(39)43-33-35(34-45-47(41,42)44-32-31-38)46-37(40)30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h13-16,35H,3-12,17-34,38H2,1-2H3,(H,41,42)"
mets_smiles["s_3427"] = "InChI=1S/C23H46NO7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-23(26)29-20-22(25)21-31-32(27,28)30-19-18-24/h9-10,22,25H,2-8,11-21,24H2,1H3,(H,27,28)/b10-9+"
mets_smiles["s_3472"] = "InChI=1S/C23H46NO7P/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-23(26)29-20-22(25)21-31-32(27,28)30-19-18-24/h9-10,22,25H,2-8,11-21,24H2,1H3,(H,27,28)/b10-9+"

for (k,v) in mets_smiles 
    if !ismissing(v) && startswith(v,"InChI")
        model.metabolites[k].annotations["InChI"] = [v]
    end
end

AA_df = DataFrame(CSV.File("data/data/yeastGEM/databases/swissprot.tsv"))

df = DataFrame(Reaction=String[],Gene=String[],enzyme=String[],substrates=String[],products=String[])
##### make a dataframe of the AA sequence, substrates and products of each reaction_gene_association
for (r,rxn) in model.reactions
    isnothing(rxn.gene_association_dnf) && continue
    for g in vcat(rxn.gene_association_dnf...)
        for row in eachrow(filter(row -> contains(row.gene_id,g),AA_df))
            subs = ";"
            prods = ";"
            for (k,v) in rxn.stoichiometry
                if haskey(model.metabolites[k].annotations,"InChI") 
                    if v < 0
                        for s in model.metabolites[k].annotations["InChI"]
                            subs = subs*s*";"
                        end
                    else
                        for p in model.metabolites[k].annotations["InChI"]
                            prods = prods*p*";"
                        end
                    end
                elseif haskey(model.metabolites[k].annotations,"smiles") 
                    if v < 0
                        for s in model.metabolites[k].annotations["smiles"]
                            subs = subs*s*";"
                        end
                    else
                        for p in model.metabolites[k].annotations["smiles"]
                            prods = prods*p*";"
                        end
                    end
                elseif haskey(model.metabolites[k].annotations,"KEGG Compound") 
                    if v < 0
                        for s in model.metabolites[k].annotations["KEGG Compound"][1]
                            subs = subs*s*";"
                        end
                    else
                        for p in model.metabolites[k].annotations["KEGG Compound"][1]
                            prods = prods*p*";"
                        end
                    end
                end
            end
            if rxn.lower_bound < 0 && rxn.upper_bound > 0 
                push!(df,(r,g,row.sequence,subs,prods))
                push!(df,("$(r)_rev",g,row.sequence,prods,subs))
            else 
                push!(df,(r,g,row.sequence,subs,prods))
            end
        end
    end
end

unique!(df)

#### save the df for using TurnUp
CSV.write("data/curated_data/yeastGEM/turnup_input.csv",df[:,[3,4,5]])
CSV.write("data/curated_data/yeastGEM/reaction_data.csv",df)

save_model(convert(JSONFBCModels.JSONFBCModel, model), "data/curated_data/yeastGEM/yeastGEMcurated.json")
