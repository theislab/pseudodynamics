function [] = mTep_fun(n)
% compute the likelihood profile for the n-th parameter using PESTO
n = str2num(n);

%% Model Definition
n_grid = 300;
x = linspace(0,1,n_grid);

%% Data
% more samples
load dataExample
clear options.u0
clear sol
load resultsExample

parameters = Est{1}.parameters;

%% Definition of the Paramter Estimation Problem

modelfun = @ simulate_pd_fv;

options.grid = x;

% compute mean u0 distribution
for i0 = 1:3
    u0(i0,:) = ksdensity(D.ind.hist{i0},x,'support',[0,1],'function','pdf');
    N01(i0) = length(D.ind.hist{i0});
    u0(i0,:) = u0(i0,:)/trapz(x,u0(i0,:));%*D.ind.size(i0);
end

options.u0 = mean(u0,1)*D.pop.mean(1);
% compute initial cells in finite volumes
options.u0 = 0.5*(options.u0(1:n_grid-1)+options.u0(2:n_grid));

options.n_exp = 1;

% compute augmentation of data ecdfs to x_combined
k1 = 1;
for it = 1:6
    options.x_combined{k1} = union(options.grid(1:end-1)',D.xsdt{it});
    D.csd_a{k1} = augment_cdf(D.xsdt{it},options.x_combined{k1},D.csdt{it});
    k1=k1+1;
end

% compute matrices for augmentation of cdfs
for it = 2:length(D.pop.t)
    for ig = 1:n_grid-1
        e_i = zeros(1,n_grid-1);
        e_i(ig) = 1;
        options.Aug_matrix{it}(:,ig) = augment_cdf(options.grid(1:end-1)',options.x_combined{it},e_i');
    end
end

% Log-likelihood function
objectiveFunction = @(theta) llPseudodynamicsFvKS(theta,modelfun,D,options);

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

save(['parametersProfileExample_' num2str(n)],'parameters')
end