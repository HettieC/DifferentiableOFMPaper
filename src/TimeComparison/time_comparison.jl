import AbstractFBCModels.CanonicalModel as CM
import DifferentiableMetabolism as D
import FastDifferentiation as F
const Ex = F.Node
import ConstraintTrees as C
import AbstractFBCModels as A
import JSONFBCModels as JFBC
import COBREXA as X
using JSON, HiGHS
using ElementaryFluxModes
using BenchmarkTools
using StatsBase
include("utils.jl")

#= Compare how long it takes to get sensitivities of all pruned reactions
to all pruned kcats using finite diff vs DiffMet =#

optimizer = HiGHS.Optimizer
gene_zero_tol = 1e-9
flux_zero_tol = 1e-9

models = [
    "data/yeastGEM/yeastGEMcurated.json",
    "data/EColi/processed_files/iML1515_processed.json"
]
isozymes = [
    "data/yeastGEM/curated/reaction_isozymes.json",
    "data/EColi/processed_files/isozymes.json"
]
masses = [
    "data/yeastGEM/curated/gene_product_molar_mass.json",
    "data/EColi/processed_files/gene_product_molar_mass.json"
]
capacities = [
    "data/yeastGEM/curated/capacity.json",
    "data/EColi/processed_files/capacity.json"
]
atps = [
    "r_4046",
    "ATPM"
]
biomasses = [
    "r_2111",
    "BIOMASS_Ec_iML1515_core_75p37M"
]
glucoses = [
    "r_1714",
    "EX_glc__D_e"
]

finitediff = Float64[]
labels = ["yeast","e. coli","yeast OFM", "e. coli OFM"]
### run finite diff on whole models
for (mod,isos,mass,cap,atp,biomass,glc) in (zip(models,isozymes,masses,capacities,atps,biomasses,glucoses))
    # run finite diff 10x 
    for i in 1:10
        fd = finite_difference(mod,isos,mass,cap,glc)
        push!(finitediff,geomean(fd.times))
    end
end
### run finite diff on ofms
finitediff_ofm = []
for (mod,isos,mass,cap,atp,biomass,glc) in (zip(models,isozymes,masses,capacities,atps,biomasses,glucoses))
    # run finite diff 10x
    for i in 1:10
        fd_ofm = finite_diff_ofms(mod,isos,mass,cap,atp,biomass,glc)
        push!(finitediff_ofm,mean(fd_ofm.times))
    end
end
all_finite_diff = vcat(finitediff,finitediff_ofm)./1E9

### get diffmet times 
diffmet_ofm = []
diffmet = []
n_reactions = Int64[]
for (mod,isos,mass,cap,atp,biomass,glc) in (zip(models,isozymes,masses,capacities,atps,biomasses,glucoses))
    for i in 1:10
        dm_efm = diff_met_ofm(mod, isos, mass, cap, atp, biomass, glc)
        append!(diffmet_ofm,dm_efm.times)
        dm,n = diff_met(mod,isos,mass,cap,glc)
        append!(diffmet,dm.times)
        append!(n_reactions,n)
    end
end
all_diff_met = vcat(diffmet,diffmet_ofm)./1E9


using DataFrames, StatsBase, CSV

## dataframe of mean and std timings 
df = DataFrame(Model=String[],FiniteDiffMean=Float64[],FiniteDiffSD=Float64[],DiffMetMean=Float64[],DiffMetSD=Float64[])
for (i,model) in enumerate(labels)
    fd_mean = geomean(all_finite_diff[10*i-9:10*i])
    fd_std = std(all_finite_diff[10*i-9:10*i])
    dm_mean = geomean(all_diff_met[10*i-9:10*i])
    dm_std = std(all_diff_met[10*i-9:10*i])
    push!(df,[model,fd_mean,fd_std,dm_mean,dm_std])
end
df.N_Parameters = vcat([unique(n_reactions),2,2]...)

CSV.write("data/results/timing.csv",df)

using CairoMakie

f = Figure(;size=(800, 600), background_color=:transparent,);
xlabs = ["GD","DiffMet","FiniteDiff"]
ax = Axis(
    f[1, 1],
    aspect=AxisAspect(1),
    ylabel="Time (s)",
    xlabelsize=35,
    ylabelsize=35,
    xticklabelsize=30,
    yticklabelsize=30,
    backgroundcolor=:transparent,
    xticks = (1:3, xlabs),
    yscale=log10,
    xgridvisible=false,
    ygridvisible=false,
    )
xs = [1,2,3]
scatter!(ax,xs,[grad[1],diffmet[1],finitediff[1]]./1E9,color=:orange,markersize=20,label="Yeast")
scatter!(ax,xs,[grad[2],diffmet[2],finitediff[2]]./1E9,color=:green,markersize=20,label="E. coli")
scatter!(ax,[2,3],[diffmet_ofm[1],finitediff_ofm[1]]./1E9,color=:orange,marker=Rect,markersize=20,label="Yeast OFM")
scatter!(ax,[2,3],[diffmet_ofm[2],finitediff_ofm[2]]./1E9,color=:green,marker=Rect,markersize=20,label="E. coli OFM")
axislegend(ax,labelsize=25,position=:lt)
f
