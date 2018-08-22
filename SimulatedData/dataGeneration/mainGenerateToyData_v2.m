clear;
n_grid=300;
x = linspace(0,1,n_grid);
u0 = 15000*lognpdf(x,-2.7,0.8);
% distribution in volumes
options.u0 = 0.5*(u0(1:n_grid-1)+u0(2:n_grid));

parameters.min = [-10.3616*ones(1,9),-11.5129*ones(1,9),-6*ones(1,9)]';
parameters.max = [0*ones(1,9),0*ones(1,9),5*ones(1,9)]';
parameters.name = {'D1','D2','D3','D4','D5','D6','D7','D8','D9',...
                'v1','v2','v3','v4','v5','v6','v7','v8','v9',...
                'a1','a2','a3','a4','a5','a6','a7','a8','a9'};
parameters.number = length(parameters.name);

% number of simulated single cells
options.nSample = 10000;
% number of replicates for single cells
options.nReplicate = 3;
% number of population size measurement replicates
options.pSample = 20;

% time points of population size measurements
D.pop.t = 0:5;
% time points of single cell measurement replicates
D.ind.tp = [0*ones(1,options.nReplicate),1*ones(1,options.nReplicate),2*ones(1,options.nReplicate),3*ones(1,options.nReplicate),4*ones(1,options.nReplicate),5*ones(1,options.nReplicate)];

%% forward simulation B
rng(2006)
% variance for data
options.varPop = [1e7,5e7,2e8,2e9,2e10,2e11];

parameters.true = [-7; -7; -7; -7; -7; -7; -7; -7; -8;... % D
    -4; -2; -2; -2; -3; -4; -4; -6; -7;... % v
    2.5;0.5;0.5;0.5;1.5;1;2;0;0]; % a

options.name = 'dataB';
generateToyData_v2(parameters.true,@simulate_pd_fv,D,options)