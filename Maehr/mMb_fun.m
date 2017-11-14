function [] = mMb_fun(n)
% compute the n-th start of a multistart optimization of the parameters for 
% the Maehr data using PESTO
n = str2num(n);

rng(10061988)

%% Model Definition
% this needs to be defined coherently with the *_syms model file
n_grid = 300;
x = linspace(0,1,n_grid);

%% Data
ic = nan;
load dataMaehr8d_b_wN

% determine crossvalidation fold and regularization alpha depending on n
if n<41
    ic = 6;
    options.alpha = 30;
end
if isnan(ic)
elseif ic == 8
    D.pop.t = D.pop.t([1:ic-1]);
    D.pop.mean = D.pop.mean([1:ic-1]);
    D.pop.var = D.pop.var([1:ic-1]);
    D.weight.b1.mean = D.weight.b1.mean([1:ic-1]);
    D.weight.b1.var = D.weight.b1.var([1:ic-1]);
    D.dist.b1.mean = D.dist.b1.mean([1:ic-1]);
    D.dist.b1.var = D.dist.b1.var([1:ic-1]);
    D.dist.b2.mean = D.dist.b2.mean([1:ic-1]);
    D.dist.b2.var = D.dist.b2.var([1:ic-1]);
else
    D.pop.t = D.pop.t([1:ic-1,ic+1:end]);
    D.pop.mean = D.pop.mean([1:ic-1,ic+1:end]);
    D.pop.var = D.pop.var([1:ic-1,ic+1:end]);
    D.weight.b1.mean = D.weight.b1.mean([1:ic-1,ic+1:end]);
    D.weight.b1.var = D.weight.b1.var([1:ic-1,ic+1:end]);
    D.dist.b1.mean = D.dist.b1.mean([1:ic-1,ic+1:end]);
    D.dist.b1.var = D.dist.b1.var([1:ic-1,ic+1:end]);
    D.dist.b2.mean = D.dist.b2.mean([1:ic-1,ic+1:end]);
    D.dist.b2.var = D.dist.b2.var([1:ic-1,ic+1:end]);
end

%% Definition of the Paramter Estimation Problem

% parameters
model = 'branching_fv';

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

% artificially shift cell at s=0 to s=e-14
D.ind.b1.hist{3}(153) = 10^(-20);

% compute mean u0 distribution
for i0 = 1:3
    u0(i0,:) = ksdensity(D.ind.b1.hist{i0},x,'support',[0,1],'function','pdf');
    N01(i0) = length(D.ind.b1.hist{i0});
    u0(i0,:) = u0(i0,:)/trapz(x,u0(i0,:));%*D.ind.size(i0);
end
for i0 = 1:3
    u02(i0,:) = ksdensity(D.ind.b2.hist{i0},x(20:end),'support',[0,1],'function','pdf');
    N02(i0) = length(D.ind.b2.hist{i0});
    u02(i0,:) = u02(i0,:)/trapz(x(20:end),u02(i0,:));%*D.ind.size(i0);
end
% weighting of subpopulations
w1=mean(N01./(N01+N02));
w2=mean(N02./(N01+N02));

options.u0 = [w1*mean(u0,1),w2*mean(u02,1)]*D.pop.mean(1);

options.u0 = 0.5*([options.u0(1:n_grid-1)+options.u0(2:n_grid),options.u0(n_grid+1:end-1)+options.u0(n_grid+2:end)]);

% number of experiments is already include in data extraction
options.n_exp = 1;

% compute augmentation of data ecdfs to x_combined
k1 = 1;
for it = 1+D.pop.t
options.x_combined{1,k1} = union(options.grid(1:end-1)',D.b1.xsdt{it});
D.b1.csd_a{k1} = augment_cdf(D.b1.xsdt{it},options.x_combined{1,k1},D.b1.csdt{it});
                           
options.x_combined{2,k1} = union(options.grid(20:end-1)',D.b2.xsdt{it});
D.b2.csd_a{k1} = augment_cdf(D.b2.xsdt{it},options.x_combined{2,k1},D.b2.csdt{it});
k1 = k1+1;
end
 
% compute matrices for augmentation of cdfs
for it = 2:length(D.pop.t)
    for ig = 1:n_grid-1
        e_i = zeros(1,n_grid-1);
        e_i(ig) = 1;
        options.Aug_matrix{1,it}(:,ig) = augment_cdf(options.grid(1:end-1)',options.x_combined{1,it},e_i');
    end
    for ig = 1:length(options.grid(20:end-1))
        e_i = zeros(1,length(options.grid(20:end-1)));
        e_i(ig) = 1;
        options.Aug_matrix{2,it}(:,ig) = augment_cdf(options.grid(20:end-1)',options.x_combined{2,it},e_i');
    end
end

                              
% Log-likelihood function
objectiveFunction = @(theta) llPseudodynamicsFvMaehrKSb4(theta,modelfun,D,options);

%% Multi-start local optimization
% A multi-start local optimization is performed within the bounds defined in
% parameters.min and .max in order to infer the unknown parameters from 
% measurement data. Therefore, a PestoOptions object is created and
% some of its properties are set accordingly.

% Options
optionsMultistart = PestoOptions();
optionsMultistart.obj_type = 'negative log-posterior';
optionsMultistart.n_starts = 120;
optionsMultistart.comp_type = 'sequential';
optionsMultistart.mode = 'silent';
optionsMultistart.proposal = 'uniform';
optionsMultistart.start_index = mod(n-1,80)+1;

% Options PESTO-application_note
optionsMultistart.fmincon.Display = 'off';
optionsMultistart.fmincon.Gradobj = 'on';
optionsMultistart.fmincon.MaxIter = 6000;
optionsMultistart.fmincon.MaxFunEvals = 12000;

% Options PESTO-master
% optionsMultistart.localOptimizerOptions.Display = 'off';
% optionsMultistart.localOptimizerOptions.Gradobj = 'on';
% optionsMultistart.localOptimizerOptions.MaxIter = 6000;
% optionsMultistart.localOptimizerOptions.MaxFunEvals = 12000;


% Optimization
parameters = getMultiStarts(parameters, objectiveFunction, optionsMultistart);

save(['parametersMaehr_' num2str(n)],'parameters','ic','options','optionsMultistart')
end