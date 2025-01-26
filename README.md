# Differentiable EFMs paper

## Description

This repository contains all required code and data to reproduce the figures in the paper 'Fast differentiation and sensitivity analysis of optimal flux modes in metabolic models'. 

### For those interested in calculating sensitivities of elementary flux modes and optimal flux modes in enzyme constrained genome scale metabolic models, please see the DifferentiableMetabolism.jl and ElementaryFluxModes.jl packages.

To reproduce the table and figures used in the paper, take the following steps:

## Table 1
Raw timings are saved in `data/results/timing.csv`. 
To calculate timings on one's one setup, run:
`src/TimeComparison/time_comparison.jl`

## Fig. 1
Genome scale models, kcats, protein fasta files all need to be downloaded from https://zenodo.org/records/6438262 (Li, F., Yuan, L., Lu, H. et al. Deep learning-based kcat prediction enables improved enzyme-constrained model reconstruction. Nat Catal 5, 662–672 (2022). https://doi.org/10.1038/s41929-022-00798-z). Next, editing file paths to `data/data/fungi343species/PredcitedKcat343species`, `data/data/fungi343species/Proteinfasta`, `data/data/fungi343species/ssGEMs` as required, run the code in the following order:
1. `src/Fungi/setup/gecko_setup.jl`
2. `src/Fungi/sensitivities.jl`
3. `src/Fungi/plotting.jl`

## Fig. 2
iML1515 model was taken from http://bigg.ucsd.edu/models/iML1515, to follow the curation steps to make an enzyme-constrained model, run:
`data/code/EColi/gecko_setup.jl`
If only interested in reproducing the figure, run:
`src/EColi/efm.jl`

## Fig. 3 
As in Figure 2, running `data/code/EColi/gecko_setup.jl` is not necessary. To reproduce the figure, run:
`src/EColi/acetate.jl`
