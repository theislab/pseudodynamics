function [] = mTp_b_fun(n)
% compute the profile likelihood for the n-th parameter for the simulated 
% branching data set using PESTO
n = str2num(n);

%% Model Definition
n_grid = 300;
x = linspace(0,1,n_grid);

%% Data

% more samples
load dataBranching
clear options.u0
clear sol
load parametersBranching
clear parameters
clear optionsMultistarts
load resultsBranching

parameters = Est{1}.parameters;

%% Definition of the Paramter Estimation Problem

modelfun = @ simulate_pd_branching_fv_toy;

k1 = 1;
for it = 1+D.pop.t
D.b1.csd_a{k1} = augment_cdf(D.b1.xsdt{it},options.x_combined{1,k1},D.b1.csdt{it});
D.b2.csd_a{k1} = augment_cdf(D.b2.xsdt{it},options.x_combined{2,k1},D.b2.csdt{it});
k1 = k1+1;
end
                              
% Log-likelihood function
objectiveFunction = @(theta) llPseudodynamicsFvKSb_toy(theta,modelfun,D,options);

%% Multi-start local optimization
% A multi-start local optimization is performed within the bounds defined in
% parameters.min and .max in order to infer the unknown parameters from 
% measurement data. Therefore, a PestoOptions object is created and
% some of its properties are set accordingly.

% Options
optionsMultistart = PestoOptions();
optionsMultistart.obj_type = 'negative log-posterior';
optionsMultistart.comp_type = 'sequential';
optionsMultistart.mode = 'text';
optionsMultistart.proposal = 'uniform';
optionsMultistart.localOptimizerOptions.Display = 'off';
optionsMultistart.localOptimizerOptions.Gradobj = 'on';
optionsMultistart.localOptimizerOptions.MaxIter = 6000;
optionsMultistart.localOptimizerOptions.MaxFunEvals = 12000;
optionsMultistart.parameter_index = n;
                               
% Optimization
parameters = getParameterProfiles(parameters, objectiveFunction, optionsMultistart);

save(['parametersProfileBranching_' num2str(n)],'parameters')
end                                