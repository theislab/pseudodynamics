# Introduction
This folder contains the scripts necessary to fit the pseudodynamics model to the data presented in Klein et al., 2015 (https://www.cell.com/cell/abstract/S0092-8674(15)00500-0), and in particular to the diffusion pseudotime representation presented in Haghverdi et al., 2016 (https://www.nature.com/articles/nmeth.3971). We published a script to generate the pseudodynamics input in manuscript accompanying this repository in Supp. data 4.1. The input is also directly supplied in Supp. data 6. We reference to this input file as pseudodynamics_input_Klein15_dpt_DPTpaperProcessedData.csv.

# Work flow
## Compute cell state ordering
1. Download the processed single-cell RNA-seq data from www.nature.com/nmeth/journal/v13/n10/extref/nmeth.3971-S2.zip
2. Extraction the processed count matrix from the matlab file with loadDataFromMatlabFile.m
3. Compute the cell states (diffusion pseudotime) using scanpy in the jupyter notebook Klein15.ipynb

The result of this workflow is pseudodynamics_input_Klein15_dpt_DPTpaperProcessedData.csv 

## Pseudodynamics
1. Compile the pseudodynamics model /pseudodynamics/finiteVolume/models/pd_fv_Klein_syms.m with pd_fv_wrap.
2. Perform the multip-start parameter estimation by running mK_fun(n) for n = 1:320.
3. Obtain the estimation results with analyzeServerResultsKlein.m.
