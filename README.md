# SeedScripts
# Scripts provided for the Paper "Contextualized metabolic modeling revealed factors affecting isoflavone accumulation in soybean seeds"
By Carolina A. Contador, Ailin Liu, Ming-Sin Ng, Yee-Shan Ku, Siu H. J. Chan, and Hon-Ming Lam 

## The provided scripts require the following tools to be available:

Cobra Toolbox ([available here](https://github.com/opencobra/cobratoolbox)) and Gurobi (Academic licenses are available free of charge via the ([Gurobi academic program](https://www.gurobi.com/academia/academic-program-and-licenses/)). The solver used in the paper was Gurobi v9.5.0.  


## Matlab Function Descriptions:
All Matlab Scripts are contained in the *Matlab* folder.
Validation.m - Model against published experimental data for C08 and W05.
ATPM_optimization.m: strategy to determine ATP and NADPH for maintenaince. C08 is used as an example.
ISO_production_no_maternal.m: isoflavonoid production in C08 and W05 seeds without maternal contribution.
post_processing.m - sampling and remove thermodynamically infeasible internal cycles in the flux distributions from the sampling.
Models folder - C08 and W05 Cotyledon models


## Soybean model:
xml files were downloaded from supplementary material of ([Moreira et al (2019)](https://doi.org/10.1104/pp.19.00122)).

## Please cite:
Contador, C.A., Liu, A., Ng, M-S., Ku, Y-S., Chan, S.H.J., Lam, H-M. Contextualized metabolic modeling revealed factors affecting isoflavone accumulation in soybean seeds. Plant Cell Environ. 2024 Sep 18.
doi: 10.1111/pce.15140. Online ahead of print. 
