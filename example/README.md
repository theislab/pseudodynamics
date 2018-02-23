# Workflow

## data generation
/pseudodynamics/simulatedData/dataGeneration/mainGenerateToyData.m

## data extraction
extractDataInput.m

## pseudodynamics
1. compile /pseudodynamics/finiteVolume/models/pd_fv_syms.m using pd_fv_wrap.m
2. evaluate mTe_fun('n') for n=1:240
3. collectParameters.m
4. sortResultsExample.m
5. analyzeResults.m
