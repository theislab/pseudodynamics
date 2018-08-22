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

% number of simulated single cells
options.nSample = 10000;
% number of population size measurement replicates
options.pSample = 20;

% variance for data
options.varPop = [1e7,5e7,2e8,2e9,2e10,2e11];

% time points of population size measurements
D.pop.t = 0:5;
% time points of single cell measurement replicates
D.ind.tp = [0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5];

% attractor
parD = -6.7;
parA = 1.2;
% true parameter
parameters.true = [parD*ones(9,1);-2;-2;-2;-2;-2;-4;-4;-10;-10;parA*ones(9,1)];

options.name = 'dataAttractor';
generateToyData_v1(parameters.true,@simulate_pd_fv,D,options)

% quasi steady state
parD = -7;
parV = -2;
% true parameter
parameters.true = [parD*ones(9,1);parV*ones(9,1);1;1;1;1.4;3.0;-0.8;-5;-5.5;-5.5];

options.name = 'dataQuasiSteady';
generateToyData_v1(parameters.true,@simulate_pd_fv,D,options)

% forward simulation
% population size variance
options.varPop = [1e6,2e6,3e6,4e6,5e6,6e6];
% true parameters
parameters.true = [-7; -8; -9; -8; -7; -7; -7; -6; -6;...
    -2; -1; -1; -1; -2; -3; -4; -4; -5;...
    1;1;0.8;0.7;0.6;0.5;0.4;0.3;0.2];

options.name = 'dataToy';
generateToyData_v1(parameters.true,@simulate_pd_fv,D,options)
