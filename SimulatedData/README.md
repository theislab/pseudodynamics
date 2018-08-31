# Introduction
We performed multiple simulation studies in our manuscript to validate different aspects of the pseudodynamics model inference.  In particual, we performed arameter estimation based on 
1. a forward simulation, with parameterization 1 (A)
2. a forward simulation, with parameterization 2 (B)
3. data generated under an attractor (C1) and a quasi steady state (C2) model. 
4. data generated with cell state sampling bias (D)
5. data generated on a branched cell state model (E).

Please also refer to our Supp. methods Sec. 8 for details on all of these case studies.

# Work flow
## Data generation
1. Simulate data with ./dataGeneration/mainGenerateToyData.m

## Pseudodynamics
### Forward simulations 1,2 (A,B)
1. Compile /pseudodynamics/finiteVolume/models/pd_fv_syms.m using pd_fv_wrap
2. Run multi-start optimzation (parameter estimation) mTf_fun('n') for n=...
3. Run multi-start optimzation (likelihood profile computation) mTfp_fun('n') for n=...

### Attractor (C1)
1. Compile /pseudodynamics/finiteVolume/models/pd_fv_syms.m using pd_fv_wrap
2. Run multi-start optimzation (parameter estimation) mTa_fun('n') for n1:240
3. Run multi-start optimzation (likelihood profile computation) mTap_fun('n') for n=...

### Quasi-steady state (C2)
1. Compile /pseudodynamics/finiteVolume/models/pd_fv_syms.m using pd_fv_wrap
2. Run multi-start optimzation (parameter estimation) mTq_fun('n') for n=1:240 
3. Run multi-start optimzation (likelihood profile computation) mTqp_fun('n') for n=...

### Sampling bias (D)
1. Compile /pseudodynamics/finiteVolume/models/pd_fv_syms.m using pd_fv_wrap
2. a) Run multi-start optimzation (parameter estimation) data simulated with bias and without correction mT_wB_fun('n') for n=... (E1)
2. b) Run multi-start optimzation (parameter estimation)for data simulated without bias and without correction mT_woB_fun('n') for n=... (E2)
2. c) Run multi-start optimzation (parameter estimation)for data simulated with bias and with correction mT_woB_fun('n') for n=... (E3)
3. a) Run multi-start optimzation (likelihood profile computation) for data simulated with bias and without correction mTp_wB_fun('n') for n=... (E1)
3. b) Run multi-start optimzation (likelihood profile computation)for data simulated without bias and without correction mTp_woB_fun('n') for n=... (E2)
3. c) Run multi-start optimzation (likelihood profile computation) for data simulated with bias and with correction mTp_woB_fun('n') for n=... (E3)

### Forward simulation of branched model (E)
1. Compile /pseudodynamics/finiteVolume/models/pd_fv_syms.m using pd_fv_wrap
2. Run multi-start optimzation mT_b_fun('n') for n=...
