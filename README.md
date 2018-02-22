# Pseudodynamics
This repository contains the example code for using pseudodynamics on all application cases presented in the associated manuscript "Beyond pseudotime: Learning population dynamics from single-cell RNAseq time series" (https://www.biorxiv.org/content/early/2017/11/14/219188.figures-only).

## Installation guide
Pseudodynamics itself is not a software package that has to be installed, but consists of a set of Matlab scripts that have to be adapted and called within a Matlab session for each application case.

Pseudodynamics depends on the forward simulation of and parameter estimation on a partial-differential equation (PDE) system capabilities of AMICI (https://github.com/ICB-DCM/AMICI) and PESTO (https://github.com/ICB-DCM/PESTO). Both are third-party Matlab interfaces/toolboxes which have to be installed for pseudodynamics to work. Please follow the installation instruction on the Github repositories to install these toolboxes.

## Run time
Parameter estimation on PDEs is computationally very intensive and we therefore recommend using pseudodynamics on a remote server and not on a local machine. Pseudodynamics parameter estimation is performed as a multi-start optimization. Therefore, the total run time is linearly dependent on the number of starts selected in the optimisation. Please refer to the supplementary methods of our manuscript for further details on the parameter estimation. We included a demo data set and code in this repository, the run time of pseudodynamics on this data set is approximately X hours on X cores for X starts in the multi-start optimisation.
