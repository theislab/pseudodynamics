function [] = mK_fun(n)
% compute the n-th start of a multistart optimization of the parameters for the Klein data using PESTO

n = str2num(n);

% set random seed
rng(10061988)

%% Model Definition
% this needs to be defined coherently with the *_syms model file
n_grid = 300;
x = linspace(0,1,n_grid);

%% Data
load dataKlein

% determine crossvalidation fold and regularization alpha depending on n
ic = nan;

if n<41
options.alpha = 1;
elseif n<81
ic = 2;
options.alpha = 1;
elseif n<121
ic = 3;
options.alpha = 1;
elseif n<161
ic = 4;
options.alpha = 1;
elseif n<201
options.alpha = 10;
elseif n<241
ic = 2;
options.alpha = 10;
elseif n<281
ic = 3;
options.alpha = 10;
elseif n<321
ic = 4;
options.alpha = 10;
end
if isnan(ic)
elseif ic == 4
D.pop.t = D.pop.t([1:ic-1]);
else
D.pop.t = D.pop.t([1:ic-1,ic+1:end]);
end


%% Definition of the Paramter Estimation Problem

% parameters
model = 'finiteVolumeKlein';

switch model
    case 'finiteVolume'
        parameters.name = {'D1','D2','D3','D4','D5','D6','D7','D8','D9',...
            'v1','v2','v3','v4','v5','v6','v7','v8','v9',...
            'a1','a2','a3','a4','a5','a6','a7','a8','a9'};
        parameters.number = length(parameters.name);
        parameters.min = -10*ones(1,parameters.number);
        parameters.max = [1*ones(1,9),0*ones(1,9),4*ones(1,9)];
        modelfun = @ simulate_pd_fv;
    case 'branching_fv'
        parameters.name = {'D1','D2','D3','D4','D5','D6','D7','D8','D9','D10','D11','D12'...
            'v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','v11','v12'...
            'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','a11','a12','d12','d21'};
        parameters.number = length(parameters.name);
        parameters.min = [-10.3616*ones(1,12),-11.5129*ones(1,12),-6*ones(1,12),-10.3616,-10.3616];
        parameters.max = [0*ones(1,12),0*ones(1,12),5*ones(1,12),4.6052,4.6052];
        modelfun = @ simulate_pd_branching_fv;
    case 'finiteVolumeKlein'
        parameters.name = {'D1','D2','D3','D4','D5','D6','D7','D8','D9'...
                'v1','v2','v3','v4','v5','v6','v7','v8','v9'};
        parameters.number = length(parameters.name);
        parameters.min = [-10.3616*ones(1,9),-11.5129*ones(1,9)];
        parameters.max = [0*ones(1,9),0*ones(1,9)];
        modelfun = @ simulate_pd_fv_Klein;
end

options.grid = x;

% artificially shift cell at s=0 to s=e-20
D.ind.hist{1}(22) = 10^(-20);

% compute u0 distribution
options.u0 = ksdensity(D.ind.hist{1},x,'support',[0,1],'function','pdf');

% compute initial cells in finite volumes
options.u0 = 0.5*([options.u0(1:n_grid-1)+options.u0(2:n_grid),options.u0(n_grid+1:end-1)+options.u0(n_grid+2:end)]);

% compute augmentation of data ecdfs to x_combined
k1 = 1;
for it = setdiff(1:4,ic)
    options.x_combined{k1} = union(options.grid(1:end-1)',D.xsdt{it});
    D.csd_a{k1} = augment_cdf(D.xsdt{it},options.x_combined{k1},D.csdt{it});
                                 
    k1 = k1+1;
end

% compute matrices for augmentation of cdfs
for it = 2:length(D.pop.t)
    for ig = 1:n_grid-1
        e_i = zeros(1,n_grid-1);
        e_i(ig) = 1;
        options.Aug_matrix{it}(:,ig) = augment_cdf(options.grid(1:end-1)',options.x_combined{it},e_i');
    end
end

% number of experiments
options.n_exp = 1;

% Log-likelihood function
objectiveFunction = @(theta) lsPseudodynamicsKlein(theta,modelfun,D,options);

%% Multi-start local optimization
% A multi-start local optimization is performed within the bounds defined in
% parameters.min and .max in order to infer the unknown parameters from 
% measurement data. Therefore, a PestoOptions object is created and
% some of its properties are set accordingly.

optionsMultistart = PestoOptions();
optionsMultistart.obj_type = 'negative log-posterior';
optionsMultistart.n_starts = 40;
optionsMultistart.comp_type = 'sequential';
optionsMultistart.mode = 'silent';
optionsMultistart.start_index = mod(n-1,40)+1;
                                   
% Options for PESTO-application_note
optionsMultistart.fmincon.Display = 'off';
optionsMultistart.fmincon.Gradobj = 'on';
optionsMultistart.fmincon.MaxIter = 2000;
optionsMultistart.fmincon.MaxFunEvals = 4000;

                                   
% Options for PESTO-master
% optionsMultistart.localOptimizerOptions.Display = 'off';
% optionsMultistart.localOptimizerOptions.Gradobj = 'on';
% optionsMultistart.localOptimizerOptions.MaxIter = 2000;
% optionsMultistart.localOptimizerOptions.MaxFunEvals = 4000;


% Optimization
parameters = getMultiStarts(parameters, objectiveFunction, optionsMultistart);

save(['parameters_Klein_' num2str(n+160)],'parameters')
end