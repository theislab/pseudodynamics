clear;
%% branching area weighted cdf
R=readtable('./data/pseudodynamics_input_maehr17_dpt_branching_TCandNK.csv');
RN=readtable('./data/pseudodynamics_input_maehr17_popsize.csv');

% number of replicates for population size
n = RN{:,4};

% population information
D.pop.t = 0:7;
D.pop.mean = RN{:,2};
D.pop.var = RN{:,3}.^2./n;

ind_b = find(strcmp(R{:,1},'branching'));
% scaling to 0.9 
R{:,6} = R{:,6}*0.9/0.55;
tau_min = min(R{ind_b,6});
tau_max = max(R{ind_b,6});
% side branch
ind2a = find(strcmp(R{:,2},'NK'));
% after branching
ind2b = find(R{:,6}>=tau_min);
% branch 2
ind2 = intersect(ind2a,ind2b);
R2 = R(ind2,:);
Rbatch2 = R2{:,5};
[Rubatch2,ia2,ic2] = unique(Rbatch2);
tp2 = R2{ia2,4};
[tp2,indT2]=sort(tp2);
Rubatch2 = Rubatch2(indT2);
D_tmp2(:,1) = R2{:,6};

for iL2 = 1:length(Rubatch2)
    indE2{iL2} = find(strcmp(R2{:,5},Rubatch2(iL2))==1);
    D.ind.b2.tp(iL2) = R2{indE2{iL2}(1),4}-12.5;
    D.ind.b2.hist{iL2} = D_tmp2(indE2{iL2});
    D.ind.b2.size(iL2) = length(indE2{iL2});
end

ind = setdiff(1:length(R{:,5}),ind2);
R1 = R(ind,:);
Rbatch1 = R1{:,5};
[Rulabels1,ia1,ic1] = unique(Rbatch1);
tp1 = R1{ia1,4};
[tp1,indT1]=sort(tp1);
Rulabels1 = Rulabels1(indT1);
D_tmp1(:,1) = R1{:,6};

for iL = 1:length(Rulabels1)
    indE1{1,iL} = find(strcmp(R1{:,5},Rulabels1(iL))==1);
    D.ind.b1.tp(iL) = R1{indE1{iL}(1),4}-12.5;
    D.ind.b1.hist{iL} = D_tmp1(indE1{iL});
    D.ind.b1.size(iL) = length(indE1{iL});
end

%% area statistic
% fraction on branches at timepoints
w1 = D.ind.b1.size./(D.ind.b1.size+D.ind.b2.size);
w2 = D.ind.b2.size./(D.ind.b1.size+D.ind.b2.size);

% empirical cdf for all replicates
for i = 1: length(D.ind.b1.hist)
    [csd1{i},xsd1{i}] = ecdf(D.ind.b1.hist{i});
    csd1{i} = csd1{i}(2:end);
    xsd1{i} = xsd1{i}(2:end);
end

% empirical cdf for pooled data of branch 1
for it = 1:length(D.pop.t)
    D_t = [];
    ind_r = find(D.ind.b1.tp == D.pop.t(it));
    % pool all data for a time point
    for ir = ind_r
        D_t = [D_t;D.ind.b1.hist{ir}];
    end
    % empirical cdf for pooled data
    [csdt1{it},xsdt1{it}] = ecdf(D_t);
    xsdt1{it} = xsdt1{it}(2:end);
    D.b1.xsdt{it} = xsdt1{it};
    csdt1{it} = csdt1{it}(2:end);
    D.weight.b1.mean(it) = mean(w1(ind_r));
    D.weight.b1.var(it) = var(w1(ind_r))/length(ind_r);
    D.b1.csdt{it} = csdt1{it};
    ih=1;
    for ir = ind_r
        % augment single replicats to obtain a cdf value for all x 
        csdr1{it}(:,ih) = augment_cdf(xsd1{ir},xsdt1{it},csd1{ir});
        % compute area between pooled and single replicat cdf
        areat1{it}(ih) = sum(diff(xsdt1{it}).*abs(csdt1{it}(1:end-1)-csdr1{it}(1:end-1,ih)));
        ih=ih+1;
    end
    D.dist.b1.mean(it) = mean(areat1{it});
    D.dist.b1.var(it) = var(areat1{it})/length(ind_r);
end

% empirical cdf for all replicates
for i = 1: length(D.ind.b2.hist)
    [csd2{i},xsd2{i}] = ecdf(D.ind.b2.hist{i});
    csd2{i} = csd2{i}(2:end);
    xsd2{i} = xsd2{i}(2:end);
end

% empirical cdf for pooled data of branch 2
for it = 1:length(D.pop.t)
    D_t = [];
    ind_r = find(D.ind.b2.tp == D.pop.t(it));
    % pool all data for a time point
    for ir = ind_r
        D_t = [D_t;D.ind.b2.hist{ir}];
    end
    % empirical cdf for pooled data
    [csdt2{it},xsdt2{it}] = ecdf(D_t);
    xsdt2{it} = xsdt2{it}(2:end);
    D.b2.xsdt{it} = xsdt2{it};
    csdt2{it} = csdt2{it}(2:end);
    D.b2.csdt{it} = csdt2{it};
    ih=1;
    for ir = ind_r
        % augment single replicats to obtain a cdf value for all x 
        csdr2{it}(:,ih) = augment_cdf(xsd2{ir},xsdt2{it},csd2{ir});
        % compute area between pooled and single replicat cdf
        areat2{it}(ih) = sum(diff(xsdt2{it}).*abs(csdt2{it}(1:end-1)-csdr2{it}(1:end-1,ih)));
        ih=ih+1;
    end
    D.dist.b2.mean(it) = mean(areat2{it});
    D.dist.b2.var(it) = var(areat2{it})/length(ind_r);
end

save('dataMaehr8d_b_wN','D');

%% monocle area weighted cdf 
R=readtable('./data/pseudodynamics_input_maehr17_monocle2_nobranching_TC.csv');
RN=readtable('./data/pseudodynamics_input_maehr17_popsize.csv');

% number of replicates for population size
n = RN{:,4};

% population information
D.pop.t = 0:7;
D.pop.mean = RN{:,2};
D.pop.var = RN{:,3}.^2./n;

% scaling to 0.9 
R{:,5} = R{:,5}*0.9/max(R{:,5});

Rbatch = R{:,4};
[Rulabels,ia,ic] = unique(Rbatch);
tp = R{ia,3};
[tp,indT]=sort(tp);
Rulabels = Rulabels(indT);
D_tmp(:,1) = R{:,5};

for iL = 1:length(Rulabels)
    indE{iL} = find(strcmp(R{:,4},Rulabels(iL))==1);
    D.ind.tp(iL) = R{indE{iL}(1),3}-12.5;
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

save('dataMaehr8d_monocle','D');
