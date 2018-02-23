# work flow

## Data generation
./dataGeneration/mainGenerateToyData.m

## Pseudodynamics
1. Compile /pseudodynamics/finiteVolume/models/pd_fv_syms.m using pd_fv_wrap
2. a) run mTa_fun('n') for n=1:240 for attractor
2. b) run mTq_fun('n') for n=1:240 for quasi steady state
