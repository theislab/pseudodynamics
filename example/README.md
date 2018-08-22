# work flow

## Data generation
/pseudodynamics/simulatedData/dataGeneration/mainGenerateToyData.m

## Data extraction
extractDataInput.m

## Pseudodynamics
1. compile /pseudodynamics/finiteVolume/models/pd_fv_syms.m using pd_fv_wrap.m
2. evaluate mTe_fun('n') for n=1:240
3. collectParameters.m
4. sortResultsExample.m
5. mTep_fun('n') for n=1:27 (profile calculation)
6. collectProfiles.m
7. analyzeResults.m
