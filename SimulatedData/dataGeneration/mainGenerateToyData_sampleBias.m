clear;
n_grid=300;
x = linspace(0,1,n_grid);
u0 = 15000*lognpdf(x,-2.7,0.8);
% distribution in volumes
options.u0 = 0.5*(u0(1:n_grid-1)+u0(2:n_grid));

parameters.name = {'D1','D2','D3','D4','D5','D6','D7','D8','D9',...
                'v1','v2','v3','v4','v5','v6','v7','v8','v9',...
                'a1','a2','a3','a4','a5','a6','a7','a8','a9'};
parameters.number = length(parameters.name);

rng(70)

% number of simulated single cells
options.nSample = 10000;
% number of replicates for single cells
options.nReplicate = 3;
% number of population size measurement replicates
options.pSample = 20;

% with bias
% options.a = 1;
% options.b = 0;

% without bias
options.a = 0;
options.b = 1;

% variance for data
options.varPop = [3e6,7e6,1e7,7e6,7e6,5e6];

% time points of population size measurements
D.pop.t = 0:5;
% time points of single cell measurement replicates
D.ind.tp = [0*ones(1,options.nReplicate),1*ones(1,options.nReplicate),2*ones(1,options.nReplicate),3*ones(1,options.nReplicate),4*ones(1,options.nReplicate),5*ones(1,options.nReplicate)];

parameters.true = [-7;-7;-7;-7;-6.5;-6.5;-6.5;-7;-8;-3;-2;-1.75;-2;-2.5;-2;-4;-5;-6;0.5;1;0.6;-0.2;-0.4;0;0.2;0.3;0.2];

options.name = 'dataSampling_woBias';
generateToyData_samplingBias(parameters.true,@simulate_pd_fv,D,options)