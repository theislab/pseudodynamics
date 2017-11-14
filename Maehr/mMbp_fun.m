function [] = mMbp_fun(n)
% compute the n-th parameter profile for the Maehr data using PESTO
n = str2num(n);

%% Model Definition
n_grid = 300;
x = linspace(0,1,n_grid);

%% Data
load dataMaehr8d_b_wN

%% load estimation result
load par10
% load par100
% load par30

%% Definition of the Paramter Estimation Problem
parameters = par{1}.parameters;
options = par{1}.options;

modelfun = @ simulate_pd_branching_fv;

% artificially shift cell at s=0 to s=e-14
D.ind.b1.hist{3}(153) = 10^(-20);

k1 = 1;
for it = 1+D.pop.t
D.b1.csd_a{k1} = augment_cdf(D.b1.xsdt{it},options.x_combined{1,k1},D.b1.csdt{it});

D.b2.csd_a{k1} = augment_cdf(D.b2.xsdt{it},options.x_combined{2,k1},D.b2.csdt{it});
k1 = k1+1;
end

% Log-likelihood function
objectiveFunction = @(theta) llPseudodynamicsFvMaehrKSb(theta,modelfun,D,options);


%% Profile likelihood calculation -- Parameters
% The uncertainty of the estimated parameters is visualized by computing
% and plotting profile likelihoods. In getParameterProfiles, this is done
% by using repeated reoptimization
% Options PESTO-master
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
                                 
parameters = getParameterProfiles(parameters, objectiveFunction, optionsMultistart);

save(['parametersProfile10_' num2str(n)],'parameters')
end