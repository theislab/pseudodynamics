# work flow
## Obtain data
1. download processed data set from www.nature.com/nmeth/journal/v13/n10/extref/nmeth.3971-S2.zip

## Compute cell state ordering
1. run data extraction script loadDataFromMatlabFile.m
2. compute diffusion pseudotime with scanpy Klein15.ipynb
3. re-format scanpy output to conform with pseudodynamics requirements createPseudodynamicsInput_Klein15.R 

## Pseudodynamics
1. pd_fv_wrap for /pseudodynamics/finiteVolume/models/pd_fv_Klein_syms.m
2. evaluate mK_fun(n) for n = 1:320
3. analyzeServerResultsKlein.m
