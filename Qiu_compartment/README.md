# Introduction
This folder contains the scripts necessary to fit the pseudodynamics model to the pancreatic beta-cell maturation data presented in our manuscript. Note that this is discussed in the context of the single-cell RNA-seq data of the same system by Qiu et al., so that we termed this study Qiu_compartment. However, no data from Qiu et al. are used here. We fit three pseudodynamics model variations to this data: a model with state-dependent transition and birth-death rates (A), a model with state-dependent transition and state- and time-dependent birth-death rates (B) and a model with state- and time-dependent transition and birth-death rates (C). Model A is nested in B and A,B are nested in C. The input is directly supplied in Supp. data 6. We reference to this input file as pseudodynamics_input_qiu17.csv and to the population size input as pseudodynamics_input_qiu17_populationsize.csv.

# Work flow
## Data extraction
Please refer to the files in our manuscript for details on the computation of cell state space. Here, we assume that processed pseudodynamics input is given and show how it is loaded in Matlab.
1. extractDataQiu_ODE.m

## Compartment model
### Constant rates
1. Compile model pseudodynamics/compartmentODE/models/compartment_syms.m using compartment_wrap.m
2. Run multi-start optimization mQc_script.m

### Constant transition and time-dependent growth rates
1. Compile model pseudodynamics/compartmentODE/models/compartment_t_syms.m using compartment_wrap.m
2. Run multi-start optimization mQct_script.m

### Time-dependent rates
1. Compile model pseudodynamics/compartmentODE/models/compartment_t4_syms.m using compartment_wrap.m
2. Run multi-start optimization mQct4_script.m
