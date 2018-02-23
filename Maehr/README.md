# work flow

## Data extraction
1. extractDataInput.m
2. extractDataRag.m

## Pseudodynamics
### Maehr Monocle
1. compile /pseudodynamics/finiteVolume/models/pd_fv_syms.m using pd_fv_wrap
2. run mMbM_fun('n') for n=1:1640

### Maehr branching
1. compile /pseudodynamics/finiteVolume/models/pd_branching_fv_syms.m using pd_fv_wrap
2. mMb_fun('n')
3. profiles using mMbp_fun.m

### Maehr Rag KO branching
1. mMbRag_fun('n')

## forward analysis for Rag KO
1. forwardAnalysisRag.m
