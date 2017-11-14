%% branching area weighted cdf - WT data only
R=readtable('./data/pseudodynamics_input_maehr17wrag_dpt_branching_TCandNK.csv');

% number of replicates for population size
n = [1,2];

% population information
D.pop.t = [2,4];

ind_R2 = find(strcmp(R{:,7},'Rag2KO'));
ind_R1 = find(strcmp(R{:,7},'Rag1KO'));
ind_R = union(ind_R1,ind_R2);
% scaling to 0.9 
R{:,6} = R{:,6}*0.9/0.48698223;
% side branch
ind_NK = find(strcmp(R{:,2},'NK'));
ind_R = setdiff(ind_R,ind_NK);

R_Rag = R(ind_R,:);
Rbatch_Rag = R_Rag{:,5};
[Rubatch_Rag,ia_R,ic_R] = unique(Rbatch_Rag);
tp = R_Rag{ia_R,4};
[tp,indT]=sort(tp,'ascend');
Rubatch_R = Rubatch_Rag(indT);
D_tmp(:,1) = R_Rag{:,6};

for iL = 1:length(Rubatch_R)
    ind{iL} = find(strcmp(R_Rag{:,5},Rubatch_R(iL))==1);
    D.ind.Rag.tp(iL) = R_Rag{ind{iL}(1),4}-12.5;
    D.ind.Rag.hist{iL} = D_tmp(ind{iL});
    D.ind.Rag.size(iL) = length(ind{iL});
end

%% area statistic

% empirical cdf for all replicates
for i = 1: length(D.ind.Rag.hist)
    [csd{i},xsd{i}] = ecdf(D.ind.Rag.hist{i});
    D.Rag.csd{i} = csd{i}(2:end);
    D.Rag.xsd{i} = xsd{i}(2:end);
end

D.Rag.csdr{1} = D.Rag.csd{1};
D.Rag.xsdr{1} = D.Rag.xsd{1};

for it = 2
    D_t = [];
    ind_r = find(D.ind.Rag.tp == D.pop.t(it));
     ih=1;
     xsdr = union(D.Rag.xsd{ind_r(1)},D.Rag.xsd{ind_r(2)});
     for ir = ind_r
        % augment single replicates to obtain a cdf value for all x 
        csdr(:,ih) = augment_cdf(D.Rag.xsd{ir},xsdr,D.Rag.csd{ir});
        ih=ih+1;
     end
     D.Rag.csdr{it} = mean(csdr');
     D.Rag.xsdr{it} = xsdr;
end

save('dataMaehr_Rag_wN','D');

%% branching area weighted cdf
clear

R=readtable('./data/pseudodynamics_input_maehr17wrag_dpt_branching_TCandNK.csv');
RN=readtable('./data/pseudodynamics_input_maehr17_popsize.csv');

% number of replicates for population size
n = RN{:,4};

% population information
D.pop.t = 0:7;
D.pop.mean = RN{:,2};
D.pop.var = RN{:,3}.^2./n;

ind_R2 = find(strcmp(R{:,7},'Rag2KO'));
ind_R1 = find(strcmp(R{:,7},'Rag1KO'));
ind_R = setdiff(1:length(R{:,7}),union(ind_R1,ind_R2));

R = R(ind_R,:);
ind_b = find(strcmp(R{:,1},'branching'));
% scaling to 0.9 
R{:,6} = R{:,6}*0.9/0.48698223;
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

save('dataMaehr8d_Rag_b_wN','D');
