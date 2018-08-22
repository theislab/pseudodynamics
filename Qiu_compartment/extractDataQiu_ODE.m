%% area weighted cdf
R=readtable('./rawData_ODE/pseudodynamics_input_lickert18_distribution.csv');
RN=readtable('./rawData_ODE/pseudodynamics_input_qiu17_populationsize.csv','TreatAsEmpty','NA');

% number of replicates for population size
n = RN{:,4};

% population information
D.pop.t = RN{:,1};

% log population
% 1st data point for initializing linear, 2nd-9th log
D.pop.mean = RN{:,6};
D.pop.var = RN{:,5}.^2;

tp = R{:,2};
D.ind.tp = tp;
D.ind.w2 = R{:,3:5};
D.mature.mean = mean(D.ind.w2,2);
D.mature.var = var(D.ind.w2,0,2)/size(D.ind.w2,2);

save('dataQiu_ODE','D');