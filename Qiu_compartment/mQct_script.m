% compute 100 starts of a multistart optimization of the parameters
% for the compartment model with constant transition and time-dependent growth
% rates using PESTO

rng(10061988)

load dataQiu_ODE
%% Definition of the Paramter Estimation Problem

% parameters
model = 'compartment_t';

switch model
    case 'finiteVolume'
            parameters.name = {'D1','D2','D3','D4','D5','D6','D7','D8','D9',...
                'v1','v2','v3','v4','v5','v6','v7','v8','v9',...
                'a1','a2','a3','a4','a5','a6','a7','a8','a9'};
            parameters.number = length(parameters.name);
            parameters.min = [-10.3616*ones(1,9),-11.5129*ones(1,9),-3*ones(1,9)];
            parameters.max = [0*ones(1,9),0*ones(1,9),1*ones(1,9)];
            modelfun = @ simulate_pd_fv;
    case 'branching_fv'
         parameters.name = {'D1','D2','D3','D4','D5','D6','D7','D8','D9','D10','D11','D12'...
            'v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','v11','v12'...
            'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','a11','a12','d12','d21'};
        parameters.number = length(parameters.name);
        parameters.min = [-10.3616*ones(1,12),-11.5129*ones(1,12),-6*ones(1,12),-10.3616,-10.3616];
        parameters.max = [0*ones(1,12),0*ones(1,12),5*ones(1,12),4.6052,4.6052];
        modelfun = @ simulate_pd_branching_fv;
    case 'compartment'
        parameters.name = {'t1' 't2' 'a1' 'a2'};
        parameters.number = length(parameters.name);
        parameters.min = [-10,-10,-3,-3];
        parameters.max = [1,1,2,2];
        modelfun = @ simulate_compartment;
    case 'compartment_t'
        parameters.name = {'t1' 't2' 'a11' 'a12' 'a13' 'a21' 'a22' 'a23'};
        parameters.number = length(parameters.name);
        parameters.min = [-10,-10,-3,-3,-3,-3,-3,-3];
        parameters.max = [1,1,2,2,2,2,2,2];
        modelfun = @ simulate_compartment_t;
end

% compute initial condition
options.u0 = exp(D.pop.mean(1))*[1-D.mature.mean(1),D.mature.mean(1)];

% regularization hyperparameter
options.alpha = 0;
% options.alpha = 1;
% options.alpha = 10;

% Log-likelihood function
objectiveFunction = @(theta) llCompartmentQiut(theta,modelfun,D,options);

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
optionsMultistart.mode = 'text';
optionsMultistart.proposal = 'uniform';
optionsMultistart.localOptimizerOptions.Display = 'off';
optionsMultistart.localOptimizerOptions.Gradobj = 'on';
optionsMultistart.localOptimizerOptions.MaxIter = 6000;
optionsMultistart.localOptimizerOptions.MaxFunEvals = 12000;
                               
% Optimization
parameters = getMultiStarts(parameters, objectiveFunction, optionsMultistart);

save(['parametersCompQiut_a' num2str(options.alpha)],'parameters','options','optionsMultistart')

%% evaluate fit
% options_sim.sensi = 1;
% options_sim.maxsteps = 1e6;
% options_sim.sx0 = zeros(2,parameters.number);
% 
% t = unique([D.pop.t;D.ind.tp]);
% [status,~,~,u,~,us] = modelfun(t,parameters.MS.par(:,2),options.u0,[],options_sim);
% 
% figure; plot(D.pop.t,exp(D.pop.mean),'-o'); hold on; plot(t,u(:,3),'-x');
% ylabel('pop size')
% xlabel('time')
% legend('data','simulation')
% figure; plot(D.ind.t,D.mature.mean,'-o'); hold on; plot(t,bsxfun(@rdivide,u(:,2),u(:,3)),'-x');
% xlabel('time')
% ylabel('fraction mature')
% legend('data','simulation')