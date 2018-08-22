function [] = mQ_s_logN_fun(n)
% compute the n-th start of a multistart optimization of the parameters
% for the model with state-dependent growth rate for the Qiu data using PESTO

n = str2num(n);

rng(10061988)

%% Model Definition
% this needs to be defined coherently with the *_syms model file
n_grid = 300;
x = linspace(0,1,n_grid);

%% Data
ic_p = nan;
ic_i = nan;

load dataQiu_logN

% determine (crossvalidation fold and) regularization alpha depending on n
if n<101
    options.alpha = 0;
elseif n<201
    options.alpha = 0.1;
elseif n<301
    options.alpha = 1;
elseif n<401
    ic_p = 6;
    ic_i = 7;
    options.alpha = 0;
elseif n<501
    ic_p = 6;
    ic_i = 7;
    options.alpha = 0.1;
elseif n<601
    ic_p = 6;
    ic_i = 7;
    options.alpha = 1;
end

if isnan(ic_p)
elseif ic_p == 6
    D.pop.t = D.pop.t([1:ic_p-1]);
    D.pop.mean = D.pop.mean([1:ic_p-1]);
    D.pop.var = D.pop.var([1:ic_p-1]);
else
    D.pop.t = D.pop.t([1:ic_p-1,ic_p+1:end]);
    D.pop.mean = D.pop.mean([1:ic_p-1,ic_p+1:end]);
    D.pop.var = D.pop.var([1:ic_p-1,ic_p+1:end]);
end

if isnan(ic_i)
elseif ic_i == 7
    D.ind.t = D.ind.t([1:ic_i-1]);
    D.ind.tp = D.ind.tp([1:ic_i-1]);
else
    D.ind.t = D.ind.t([1:ic_i-1,ic_i+1:end]);
    D.ind.tp = D.ind.tp([1:ic_i-1,ic_i+1:end]);
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
            parameters.min = [-10.3616*ones(1,9),-11.5129*ones(1,9),-3*ones(1,9)];
            parameters.max = [0*ones(1,9),0*ones(1,9),1*ones(1,9)];
            modelfun = @ simulate_pd_fv;
    case 'finiteVolume_t'
            parameters.name = {'D1','D2','D3','D4','D5','D6','D7','D8','D9',...
                'v1','v2','v3','v4','v5','v6','v7','v8','v9',...
                'a1','a2','a3','a4','a5','a6','a7','a8','a9'};
            parameters.number = length(parameters.name);
            parameters.min = [-10.3616*ones(1,9),-11.5129*ones(1,9),-3*ones(1,9)];
            parameters.max = [0*ones(1,9),0*ones(1,9),1*ones(1,9)];
            modelfun = @ simulate_pd_fv_t;
    case 'branching_fv'
         parameters.name = {'D1','D2','D3','D4','D5','D6','D7','D8','D9','D10','D11','D12'...
            'v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','v11','v12'...
            'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','a11','a12','d12','d21'};
        parameters.number = length(parameters.name);
        parameters.min = [-10.3616*ones(1,12),-11.5129*ones(1,12),-6*ones(1,12),-10.3616,-10.3616];
        parameters.max = [0*ones(1,12),0*ones(1,12),5*ones(1,12),4.6052,4.6052];
        modelfun = @ simulate_pd_branching_fv2;
end

options.grid = x;

% artificially shift cell at s=0 to s=e-14
D.ind.hist{1}(58) = 10^(-20);

% compute u0 distribution
options.u0 = ksdensity(D.ind.hist{1},x,'support',[0,1],'function','pdf');
options.u0 = options.u0/trapz(x,options.u0)*exp(D.pop.mean(1));

% compute initial cells in finite volumes
options.u0 = 0.5*(options.u0(1:n_grid-1)+options.u0(2:n_grid));

% compute augmentation of data ecdfs to x_combined
k1 = 1;
for it = setdiff(1:7,ic_i)
    options.x_combined{k1} = union(options.grid(1:end-1)',D.xsdt{it});
    D.csd_a{k1} = augment_cdf(D.xsdt{it},options.x_combined{k1},D.csdt{it});                                 
    k1 = k1+1;
end

% compute matrices for augmentation of cdfs
for it = 2:length(D.ind.tp)
    for ig = 1:n_grid-1
        e_i = zeros(1,n_grid-1);
        e_i(ig) = 1;
        options.Aug_matrix{it}(:,ig) = augment_cdf(options.grid(1:end-1)',options.x_combined{it},e_i');
    end
end


% Log-likelihood function
objectiveFunction = @(theta) lsPseudodynamicsQiuLogN(theta,modelfun,D,options);

%% Multi-start local optimization
% A multi-start local optimization is performed within the bounds defined in
% parameters.min and .max in order to infer the unknown parameters from 
% measurement data. Therefore, a PestoOptions object is created and
% some of its properties are set accordingly.

% Options
optionsMultistart = PestoOptions();
optionsMultistart.obj_type = 'negative log-posterior';
optionsMultistart.n_starts = 100;
optionsMultistart.comp_type = 'sequential';
optionsMultistart.mode = 'silent';
optionsMultistart.proposal = 'uniform';
optionsMultistart.localOptimizerOptions.algorithm = 'trust-region-reflective';
optionsMultistart.localOptimizerOptions.Display = 'off';
optionsMultistart.localOptimizerOptions.Gradobj = 'on';
optionsMultistart.localOptimizerOptions.MaxIter = 6000;
optionsMultistart.localOptimizerOptions.MaxFunEvals = 12000;
optionsMultistart.start_index = mod(n-1,100)+1;
                              
% Optimization
parameters = getMultiStarts(parameters, objectiveFunction, optionsMultistart);

if n==1
    save(['parametersQiu_logN_' num2str(n)],'parameters','options','optionsMultistart')
elseif n==101
    save(['parametersQiu_logN_' num2str(n)],'parameters','options','optionsMultistart')
elseif n==201
    save(['parametersQiu_logN_' num2str(n)],'parameters','options','optionsMultistart')
elseif n==301
    save(['parametersQiu_logN_' num2str(n)],'parameters','options','optionsMultistart')
elseif n==401
    save(['parametersQiu_logN_' num2str(n)],'parameters','options','optionsMultistart')
elseif n==501
    save(['parametersQiu_logN_' num2str(n)],'parameters','options','optionsMultistart')
else
    save(['parametersQiu_logN_' num2str(n)],'parameters')
end

end