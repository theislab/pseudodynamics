# work flow

## Data extraction
extractDataQiu.m

## Pseudodynamics
### Stated-dependent grwoth
1. compile pseudodynamics/finiteVolume/models/pd_fv_syms.m using pd_fv_wrap.m
2. run mQ_s_logN_fun('n') for n=1:300

### Time-dependent growth
1. compile pseudodynamics/finiteVolume/models/pd_fv_st_syms.m using pd_fv_wrap.m
2. run mQ_st_logN_fun('n') for n=1:300
