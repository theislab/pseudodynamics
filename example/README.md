# Introduction
This is a small toy example case which we supply for interested users to familiarise themselves with how pseudodynamics works. Please also refer to our vignette in `~/PseudodynamicsVignette.pdf`.

# Work flow
## Data generation
1. Simulate data to fit model to /pseudodynamics/simulatedData/dataGeneration/mainGenerateToyData.m

## Data extraction
1. Format data to comply with requirements of fitting scripts extractDataInput.m

## Pseudodynamics
1. Compile model /pseudodynamics/finiteVolume/models/pd_fv_syms.m using pd_fv_wrap.m
2. Run multi-start optimization (parameter estimation) mTe_fun('n') for n=1:240
3. Collect results of parameter estimation collectParameters.m
4. Summarize parameter estimation runs sortResultsExample.m
5. Run multi-start optimization (likelihood profile computation) mTep_fun('n') for n=1:27 (profile calculation)
6. Collect results of likelihood profile computation collectProfiles.m
7. Summarize all results analyzeResults.m
