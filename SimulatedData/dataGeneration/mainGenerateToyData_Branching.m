clear;
n_grid1=300;
n_grid2=231;
x1 = linspace(0,1,n_grid1);
x2 = linspace(0,(n_grid1-1)/(n_grid2-1),n_grid2);
u0_1 = 15000*lognpdf(x1,-2.7,0.8);
u0_2 = 2000*lognpdf(x2,-2.4,0.9);

rng(2505)
% distribution in volumes
options.u0 = [0.5*(u0_1(1:n_grid1-1)+u0_1(2:n_grid1)),0.5*(u0_2(1:n_grid2-1)+u0_2(2:n_grid2))];

parameters.name = {'D1','D2','D3','D4','D5','D6','D7','D8','D9','D10','D11','D12',...
                'v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','v11','v12',...
                'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','a11','a12','delta1_2','delta2_1'};
parameters.number = length(parameters.name);

% number of simulated single cells
options.nSample = 10000;
% number of replicates for single cells
options.nReplicate = 3;
% number of population size measurement replicates
options.pSample = 20;

% variance for data
options.varPop = [1e7,5e7,2e8,2e9,2e10,2e11];

% time points of population size measurements
D.pop.t = 0:5;
% time points of single cell measurement replicates
D.ind.b1.tp = [0*ones(1,options.nReplicate),1*ones(1,options.nReplicate),2*ones(1,options.nReplicate),3*ones(1,options.nReplicate),4*ones(1,options.nReplicate),5*ones(1,options.nReplicate)];
D.ind.b2.tp = [0*ones(1,options.nReplicate),1*ones(1,options.nReplicate),2*ones(1,options.nReplicate),3*ones(1,options.nReplicate),4*ones(1,options.nReplicate),5*ones(1,options.nReplicate)];

parameters.true = [-7;-7;-7;-7;-7;-7;-7;-7;-7;... %Db1
    -7;-7;-7;... %Db2
    -2.3;-2.3;-2.3;-2.3;-2.3;-2.3;-3;-4;-5;... %vb1
    -2.3;-3;-5;... %vb2
    1;1.2;1;0.5;0.5;0.5;0.5;-1.5;-2;... %gb1
    0.8;1.2;1;... %gb2
    0;-7];

options.name = 'dataBranching';
generateToyData_branching(parameters.true,@simulate_pd_branching_fv_toy,D,options)