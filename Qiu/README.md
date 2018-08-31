# Introduction
This folder contains the scripts necessary to fit the pseudodynamics model to the data presented in Qiu et al., 2017 (https://www.cell.com/cell-metabolism/abstract/S1550-4131(17)30208-5). We fit two pseudodynamics model variations to this data: a model with state-dependent birth-death rates only (A) and one with state- and time-dependent birth-death rates (B). We published a script to generate the pseudodynamics input in manuscript accompanying this repository in Supp. data 4.4. The input is also directly supplied in Supp. data 6. We reference to this input file as pseudodynamics_input_qiu17.csv and to the population size input as pseudodynamics_input_qiu17_populationsize.csv.

# Work flow
## Data extraction
Please refer to the files in our manuscript for details on the computation of cell state space. Here, we assume that processed pseudodynamics input is given and show how it is loaded in Matlab.
1. extractDataQiu.m

## Pseudodynamics
### Stated-dependent birth-death rates (A)
1. Compile model pseudodynamics/finiteVolume/models/pd_fv_syms.m using pd_fv_wrap.m
2. Run multi-start optimization mQ_s_logN_fun('n') for n=1:300

### State- and time-dependent birth-death rates (B)
1. Compile model pseudodynamics/finiteVolume/models/pd_fv_st_syms.m using pd_fv_wrap.m
2. Run multi-start optimization mQ_st_logN_fun('n') for n=1:300
