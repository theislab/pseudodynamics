# Introduction
This folder contains the scripts necessary to fit the pseudodynamics model to the T-cell maturation data presented in Kernfeld et al., 2018, and Fischer et al, 2018. In the manuscript, we fit the pseudodynamics model to the wild type data represented by diffusion pseudotime (A), monocle2 pseudotime (B) and to the union of wild-type, Rag1-KO and Rag2-KO data (C). A,C are branching pseudodynamics models with two branches, B only contains a single branch. A,B,C represent separate models fits. In addition, we performed a mapping of the position of beta-selection onto the T-cell maturation cell state trajectory via a separate optimization (C2). We published a script to generate the pseudodynamics input in manuscript accompanying this repository in Supp. data 4.2 and 4.3. The input is also directly supplied in Supp. data 6. 

# Work flow
## Data extraction
Please refer to the files in our manuscript for details on the computation of cell state space. Here, we assume that processed pseudodynamics input is given and show how it is loaded in Matlab.
1. extractDataInput.m
2. extractDataRag.m

## Pseudodynamics
### Diffusion pseudotime cell state space, branched pseudodynamics model, wild-type data (A)
1. Compile model /pseudodynamics/finiteVolume/models/pd_branching_fv_syms.m using pd_fv_wrap
2. Run multi-start optimization mMb_fun('n')
3. Compute parameter uncertainties via liklihood profiles using mMbp_fun.m

### Monocle2 cell state space, linear pseudodynamics model, wild-type data (B)
1. compile /pseudodynamics/finiteVolume/models/pd_fv_syms.m using pd_fv_wrap
2. Run multi-start optimization run mMbM_fun('n') for n=1:1640

### Diffusion pseudotime cell state space, branched pseudodynamics model, wild-type and knock-out data (C)
1. Run multi-start optimization mMbRag_fun('n')

#### Mapping of beta-selection check-point onto trajectory (C2)
1. Compute likelihood profile of beta-selection position forwardAnalysisRag.m
