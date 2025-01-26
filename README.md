# Differentiable EFMs paper

## Description

This repository contains all required code and data to reproduce the figures in the paper 'Fast differentiation and sensitivity analysis of optimal flux modes in metabolic models'. 

### For those interested in calculating sensitivities of elementary flux modes and optimal flux modes in enzyme constrained genome scale metabolic models, please see the DifferentiableMetabolism.jl and ElementaryFluxModes.jl packages.

To reproduce the table and figures used in the paper, follow the following steps:

## Table 1
TODO

## Fig. 1
Order of files to run:
1. `src/Fungi/getting_data/gecko_setup.jl`
2. `src/Fungi/getting_data/sensitivities.jl`
3. `src/Fungi/plotting.jl`

## Fig. 2
Order of files to run:
1. `src/EColi/gecko_setup.jl`
2. `src/EColi/efm.jl`

## Fig. 3 
Order of files to run (1. and 2. not required to re-run if Fig. 2 already produced)
1. `src/EColi/gecko_setup.jl`
2. `src/EColi/efm.jl`
3. `src/EColi/acetate.jl`