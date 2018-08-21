%% area weighted cdf
R=readtable('./rawData_new/pseudodynamics_input_qiu17.csv');
RN=readtable('./rawData_new/pseudodynamics_input_qiu17_populationsize.csv','TreatAsEmpty','NA');

% number of replicates for population size
n = RN{:,4};

% population information
D.pop.t = RN{:,1};

% log population
% 1st data point for initializing linear, 2nd-9th log
D.pop.mean = RN{:,6};
D.pop.var = RN{:,5}.^2;

% scaling to 0.9 
R{:,4} = R{:,4}*0.9/max(R{:,4});

Rbatch = R{:,3};
[Rubatch,ia,ic] = unique(Rbatch);
tp = R{ia,2};
[tp,indT]=sort(tp,'ascend');
Rubatch = Rubatch(indT);
D_tmp(:,1) = R{:,4};

for iL = 1:length(Rubatch)
    ind{iL} = find(strcmp(R{:,3},Rubatch(iL))==1);
    D.ind.tp(iL) = R{ind{iL}(1),2};
    D.ind.hist{iL} = D_tmp(ind{iL});
    D.ind.size(iL) = length(ind{iL});
end

D.ind.t = D.ind.tp';
%% area statistic

% empirical cdf for all replicates
for i = 1: length(D.ind.hist)
    [csd{i},xsd{i}] = ecdf(D.ind.hist{i});
    D.csdt{i} = csd{i}(2:end);
    D.xsdt{i} = xsd{i}(2:end);
end

save('dataQiu_logN','D');