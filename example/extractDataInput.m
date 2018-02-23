%% example simulated data
R=readtable('./rawData/pseudodynamics_input_example.csv');
RN=readtable('./rawData/pseudodynamics_input_example_popsize.csv');

% number of replicates for population size
n = RN{:,4};

% population information
D.pop.t = RN{:,1};
D.pop.mean = RN{:,2};
D.pop.var = RN{:,3}.^2./n;

% scaling to 0.9 (omitted for simulated data)
% R{:,5} = R{:,5}*0.9/max(R{:,5});

Rbatch = R{:,2};
[Rulabels,ia,ic] = unique(Rbatch);
tp = R{ia,1};
% sort experiment labels according to experiment time
[tp,indT]=sort(tp);
Rulabels = Rulabels(indT);
D_tmp(:,1) = R{:,3};

for iL = 1:length(Rulabels)
    indE{iL} = find(strcmp(R{:,2},Rulabels(iL))==1);
    D.ind.tp(iL) = R{indE{iL}(1),1};
    D.ind.hist{iL} = D_tmp(indE{iL});
    D.ind.size(iL) = length(indE{iL});
end

%% area statistic

% empirical cdf for all replicates
for i = 1: length(D.ind.hist)
    [csd{i},xsd{i}] = ecdf(D.ind.hist{i});
    csd{i} = csd{i}(2:end);
    xsd{i} = xsd{i}(2:end);
end

% empirical cdf for pooled data of branch 1
for it = 1:length(D.pop.t)
    D_t = [];
    ind_r = find(D.ind.tp == D.pop.t(it));
    % pool all data for a time point
    for ir = ind_r
        D_t = [D_t;D.ind.hist{ir}];
    end
    % empirical cdf for pooled data
    [csdt{it},xsdt{it}] = ecdf(D_t);
    xsdt{it} = xsdt{it}(2:end);
    D.xsdt{it} = xsdt{it};
    csdt{it} = csdt{it}(2:end);
    % weighted cdf
    D.csdt{it} = csdt{it};
    ih=1;
    for ir = ind_r
        % augment single replicats to obtain a cdf value for all x 
        csdr{it}(:,ih) = augment_cdf(xsd{ir},xsdt{it},csd{ir});
        % compute area between pooled and single replicat cdf
        areat{it}(ih) = sum(diff(xsdt{it}).*abs(csdt{it}(1:end-1)-csdr{it}(1:end-1,ih)));
        ih=ih+1;
    end
    D.dist.mean(it) = mean(areat{it});
    D.dist.var(it) = var(areat{it})/length(ind_r);
end

save('dataExample','D');
