# Pseudodynamics
This repository contains the example code for using pseudodynamics on all application cases presented in the associated manuscript "Beyond pseudotime: Learning population dynamics from single-cell RNAseq time series" (https://www.biorxiv.org/content/early/2017/11/14/219188.figures-only).

## Installation guide
Pseudodynamics itself is not a software package that has to be installed, but consists of a set of Matlab scripts that have to be adapted and called within a Matlab session for each application case.

Pseudodynamics depends on the forward simulation of and parameter estimation on a partial-differential equation (PDE) system capabilities of AMICI (https://github.com/ICB-DCM/AMICI) and PESTO (https://github.com/ICB-DCM/PESTO). Both are third-party Matlab interfaces/toolboxes which have to be installed for pseudodynamics to work. Please follow the installation instruction on the Github repositories to install these toolboxes.

## Examples
This repository contains the scripts necessary to perform the pseudodynamics model parameter estimation and associated workflows for all examples presented the manuscript "Beyond pseudotime: Learning population dynamics from single-cell RNAseq time series". Note that we added an introductory example as (6. example). The examples are:
1. Klein: Data from Klein et al., 2015, mouse embryonic stem cell differentiation in vitro (single-cell RNA-seq).
2. Maehr: Data from Kernfeld et al., 2018, and Fischer et al., 2018, thymic T-cell maturation (single-cell RNA-seq). This inlcudes the wild-type data represented by diffusion pseudotime and by monocle2 pseudotime.
3. Qiu: Data from Qiu et al., 2017, and Fischer et al., 2018, pancreatic beta-cell maturation (single-cell RNA-seq).
4. Qiu_compartment: Data from Fischer et al., 2018, pancreatic beta-cell maturation quantified by Ucn3 expression (cell type fractions from stained tissue sections). This exampled contains the two-compartment ordinary differential equation model for pancreatic beta-cell maturation, (3. Qiu) contains the continuous partial differential equation model.
5. SimulatedData: Data from Fischer et al., 2018, simulated toy data used for validation in manuscript.
6. example: Small example case study to familiarise users with the work flow.

The models associated with these examples are in 
- finiteVolume: continuous cell state models (partial differential equations, PDEs)
- compartmentODE: discrete cell state models (ordinary differential equations, ODEs)

## Run time
Parameter estimation on PDEs is computationally very intensive and we therefore recommend using pseudodynamics on a remote server and not on a local machine. Pseudodynamics parameter estimation is performed as a multi-start optimization. Therefore, the total run time is linearly dependent on the number of starts selected in the optimisation. Please refer to the supplementary methods of our manuscript for further details on the parameter estimation. We included a demo data set and code in this repository, the run time of pseudodynamics on this data set is on average approximately 11 hours per start on 1 core in the multi-start optimisation.
