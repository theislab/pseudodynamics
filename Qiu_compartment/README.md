# work flow

## Data extraction
extractDataQiu_ODE.m

## Compartment model
### Constant rates
1. compile pseudodynamics/compartmentODE/models/compartment_syms.m using compartment_wrap.m
2. run mQc_script.m

### Constant transition and time-dependent growth rates
1. compile pseudodynamics/compartmentODE/models/compartment_t_syms.m using compartment_wrap.m
2. run mQct_script.m

### Time-dependent rates
1. compile pseudodynamics/compartmentODE/models/compartment_t4_syms.m using compartment_wrap.m
2. run mQct4_script.m
