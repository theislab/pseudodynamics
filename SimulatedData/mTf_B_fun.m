function [] = mTf_B_fun(n)
% compute the n-th start of a multistart optimization of the parameters
% for the simulated data dataB using PESTO
n = str2num(n);

rng(10061988)

%% Model Definition
n_grid = 300;
x = linspace(0,1,n_grid);

%% Data
ic = nan;

% more samples
load dataB
clear options.u0
clear sol

% determine regularization alpha depending on n
if n<81
    options.alpha = 0;
elseif n<161
    options.alpha = 1;
elseif n<241
    options.alpha = 10;
end

%% Definition of the Paramter Estimation Problem

% parameters
model = 'finiteVolume';

switch model
    case 'finiteVolume'
            parameters.name = {'D1','D2','D3','D4','D5','D6','D7','D8','D9',...
                'v1','v2','v3','v4','v5','v6','v7','v8','v9',...
                'a1','a2','a3','a4','a5','a6','a7','a8','a9'};
            parameters.number = length(parameters.name);
            parameters.min = [-10.3616*ones(1,9),-11.5129*ones(1,9),-6*ones(1,9)];
            parameters.max = [0*ones(1,9),0*ones(1,9),5*ones(1,9)];
            parD = -6.7;
            parA = 1.2;
            parameters.guess = [parD*ones(9,1);-2;-2;-2;-2;-2;-4;-4;-10;-12;parA*ones(9,1)];
            modelfun = @ simulate_pd_fv;
    case 'branching_fv'
         parameters.name = {'D1','D2','D3','D4','D5','D6','D7','D8','D9','D10','D11','D12'...
            'v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','v11','v12'...
            'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','a11','a12','d12','d21'};
        parameters.number = length(parameters.name);
        parameters.min = [-10.3616*ones(1,12),-11.5129*ones(1,12),-6*ones(1,12),-10.3616,-10.3616];
        parameters.max = [0*ones(1,12),0*ones(1,12),5*ones(1,12),4.6052,4.6052];
        modelfun = @ simulate_pd_branching_fv;
end

options.grid = x;

% compute mean u0 distribution
for i0 = 1:3
    u0(i0,:) = ksdensity(D.ind.hist{i0},x,'support',[0,1],'function','pdf');
    N01(i0) = length(D.ind.hist{i0});
    u0(i0,:) = u0(i0,:)/trapz(x,u0(i0,:));%*D.ind.size(i0);
end

options.u0 = mean(u0,1)*D.pop.mean(1);
options.u0 = 0.5*(options.u0(1:n_grid-1)+options.u0(2:n_grid));

options.n_exp = 1;

% compute augmentation of data ecdfs to x_combined
k1 = 1;
for it = 1+D.pop.t
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
objectiveFunction = @(theta) llPseudodynamicsFvMaehrMonocleKS(theta,modelfun,D,options);

%% Multi-start local optimization
% A multi-start local optimization is performed within the bounds defined in
% parameters.min and .max in order to infer the unknown parameters from 
% measurement data. Therefore, a PestoOptions object is created and
% some of its properties are set accordingly.

% Options
optionsMultistart = PestoOptions();
optionsMultistart.obj_type = 'negative log-posterior';
optionsMultistart.n_starts = 80;
optionsMultistart.comp_type = 'sequential';
optionsMultistart.mode = 'silent';
optionsMultistart.proposal = 'uniform';
optionsMultistart.localOptimizerOptions.Display = 'off';
optionsMultistart.localOptimizerOptions.Gradobj = 'on';
optionsMultistart.localOptimizerOptions.MaxIter = 6000;
optionsMultistart.localOptimizerOptions.MaxFunEvals = 12000;
optionsMultistart.start_index = mod(n-1,80)+1;

                               
% Optimization
parameters = getMultiStarts(parameters, objectiveFunction, optionsMultistart);

if n==1
    save(['parametersB_' num2str(n)],'parameters','options','optionsMultistart')
elseif n==81
    save(['parametersB_' num2str(n)],'parameters','options','optionsMultistart')
elseif n==161
    save(['parametersB_' num2str(n)],'parameters','options','optionsMultistart')
else
    save(['parametersB_' num2str(n)],'parameters')
end
end
                                 