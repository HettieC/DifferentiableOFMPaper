# Differentiable EFMs paper

## Description

This repository contains all required code and data to reproduce the figures in the paper 'Fast differentiation and sensitivity analysis of optimal flux modes in metabolic models'. 

### For those interested in calculating sensitivities of elementary flux modes and optimal flux modes in enzyme constrained genome scale metabolic models, please see the DifferentiableMetabolism.jl and ElementaryFluxModes.jl packages.

To reproduce the table and figures used in the paper, take the following steps:

## Download required data

### yeast-GEM
- Download yeast-GEM v9.0.1 from `https://github.com/SysBioChalmers/yeast-GEM/archive/refs/tags/v9.0.1.zip`
- Extract the `yeast-GEM.mat` file into the folder `data/data/yeastGEM`.

### Fungal models
- Download https://zenodo.org/records/6438262 (Li, F., Yuan, L., Lu, H. et al. Deep learning-based kcat prediction enables improved enzyme-constrained model reconstruction. Nat Catal 5, 662–672 (2022))
- Only the folders `ssGEMs`, `PredcitedKcat343species`, and `Proteinfasta` are required

The file paths to this dataset must be edited to match:
`~/DifferentiableOFMPaper/data/data/fungi343species/ssGEMs` 
`~/DifferentiableOFMPaper/data/data/fungi343species/PredcitedKcat343species`, 
`~/DifferentiableOFMPaper/data/data/fungi343species/Proteinfasta`, 

## Reproducing figures

### Table 1
Raw timings calculated on an AMD Ryzen 9 5950X with 32 GB main memory are saved in `data/results/timing.csv`. 
To calculate timings on one's own setup, run the file:
`src/TimeComparison/time_comparison.jl`

### Fig. 1
Run the .jl files in the following order:
1. `src/Fungi/setup/gecko_setup.jl` : creates enzyme-constrained models of the ssGEMs
2. `src/Fungi/sensitivities.jl` : calculate the sensitivity of all model growth rates to their parameters
3. `src/Fungi/plotting.jl` : plot the growth sensitivity per metabolic pathway

### Fig. 2
iML1515 model was taken from http://bigg.ucsd.edu/models/iML1515, to follow the curation steps to make an enzyme-constrained model, run:
`data/code/EColi/gecko_setup.jl`
If only interested in reproducing the figure, run:
`src/EColi/efm.jl`

### Fig. 3 
As in Figure 2, running `data/code/EColi/gecko_setup.jl` is not necessary. To reproduce the figure, run:
`src/EColi/acetate.jl`
