function [] = generateToyData_branching(varargin)
% set model specifications
n_grid1 = 300;
n_grid2 = 231;
x = linspace(0,1,n_grid1);
x2 = linspace(69/299,1,n_grid2);

grid_x = linspace(x(1)+(x(2)-x(1))/2,x(end-1)+(x(end)-x(end-1))/2,n_grid1-1); 

% assign inputs
parameters_true = varargin{1};
modelfun = varargin{2};
D = varargin{3};
options = varargin{4};

% simulate model
sol = modelfun(D.pop.t,parameters_true,options.u0);
sol.p1 = bsxfun(@rdivide,sol.y(:,1:n_grid1-1),sol.y(:,end-1))*(1/(n_grid1-1));
sol.p2 = bsxfun(@rdivide,sol.y(:,n_grid1:end-2),sol.y(:,end))*(1/(n_grid1-1));

% generate number and fraction of cells on branch 1 by sampling from a
% binomial ditribution
for iT=1:length(D.pop.t)
    w(iT) = sol.y(iT,end-1)/(sol.y(iT,end-1)+sol.y(iT,end));
    D.ind.b1.size((iT-1)*options.nReplicate+1:iT*options.nReplicate) = ...
        random('bino',options.nSample,w(iT),[1,options.nReplicate]);
    D.ind.b1.weight((iT-1)*options.nReplicate+1:iT*options.nReplicate) =...
        D.ind.b1.size((iT-1)*options.nReplicate+1:iT*options.nReplicate)./...
        (options.nSample);
    D.weight.b1.mean(iT) = mean(D.ind.b1.weight((iT-1)*...
        options.nReplicate+1:iT*options.nReplicate));
    D.weight.b1.var(iT) = var(D.ind.b1.weight((iT-1)*...
        options.nReplicate+1:iT*options.nReplicate))/options.nReplicate;
end

% generate count data from distribution
% branch 1
cdf1_sim = cumsum(sol.p1,2);
k=1;
for iT = D.ind.b1.tp
    [cdf1,mask1] = unique(cdf1_sim(D.pop.t == iT,:));
    xm1 = x(mask1);
    % generate uniform random number
    rS1 = rand(1,D.ind.b1.size(k));
    % draw uniform number in interval
    for iS=1:D.ind.b1.size(k)
        try
           lb1 = xm1(find(cdf1>rS1(iS),1,'first'));
           D.ind.b1.hist{k}(iS) = lb1+1/299*rand(1);
        catch
           lb1 = 0; 
        end
    end
    k=k+1;
end

% branch 2
cdf2_sim = cumsum(sol.p2,2);
k=1;
for iT = D.ind.b2.tp
    [cdf2,mask2] = unique(cdf2_sim(D.pop.t == iT,:));
    xm2 = x2(mask2);
    % generate uniform random number
    rS2 = rand(1,options.nSample-D.ind.b1.size(k));
    % draw uniform number in interval
    for iS=1:options.nSample-D.ind.b1.size(k)
        lb2 = xm2(find(cdf2>rS2(iS),1,'first'));
        D.ind.b2.hist{k}(iS) = lb2+1/299*rand(1);
    end
    k=k+1;
end

% generate normal distributed population size measurements
for iT=1:length(D.pop.t)
    D.ind.size{iT} = normrnd(sol.y(iT,end-1)+sol.y(iT,end),...
        sqrt(options.varPop(iT)),[1,options.pSample]);
    D.pop.mean(iT) = mean(D.ind.size{iT});
    D.pop.var(iT) = var(D.ind.size{iT})/options.pSample;
end

%% area statistic
% empirical cdf for branch 1 for all replicates
for i = 1: length(D.ind.b1.hist)
    [csd1{i},xsd1{i}] = ecdf(D.ind.b1.hist{i});
    csd1{i} = csd1{i}(2:end);
    xsd1{i} = xsd1{i}(2:end);
end

% empirical cdf for pooled data of branch 1
for iT = 1:length(D.pop.t)
    D_t = [];
    ind_r = find(D.ind.b1.tp == D.pop.t(iT));
    % pool all data for a time point
    for ir = ind_r
        D_t = [D_t;D.ind.b1.hist{ir}(:)];
    end
    % empirical cdf for pooled data
    [csdt1{iT},xsdt1{iT}] = ecdf(D_t);
    xsdt1{iT} = xsdt1{iT}(2:end);
    D.b1.xsdt{iT} = xsdt1{iT};
    csdt1{iT} = csdt1{iT}(2:end);
    D.b1.csdt{iT} = csdt1{iT};
    ih=1;
    for ir = ind_r
        % augment single replicats to obtain a cdf value for all x 
        csdr1{iT}(:,ih) = augment_cdf(xsd1{ir},xsdt1{iT},csd1{ir});
        % compute area between pooled and single replicat cdf
        D.areat1{iT}(ih) = sum(diff(xsdt1{iT}).*...
            abs(csdt1{iT}(1:end-1)-csdr1{iT}(1:end-1,ih)));
        ih=ih+1;
    end
    D.dist.b1.mean(iT) = mean(D.areat1{iT});
    D.dist.b1.var(iT) = var(D.areat1{iT})/length(ind_r);
end

% empirical cdf for branch 2 for all replicates
for i = 1: length(D.ind.b2.hist)
    [csd2{i},xsd2{i}] = ecdf(D.ind.b2.hist{i});
    csd2{i} = csd2{i}(2:end);
    xsd2{i} = xsd2{i}(2:end);
end

% empirical cdf for pooled data of branch 2
for iT = 1:length(D.pop.t)
    D_t = [];
    ind_r = find(D.ind.b2.tp == D.pop.t(iT));
    % pool all data for a time point
    for ir = ind_r
        D_t = [D_t;D.ind.b2.hist{ir}(:)];
    end
    % empirical cdf for pooled data
    [csdt2{iT},xsdt2{iT}] = ecdf(D_t);
    xsdt2{iT} = xsdt2{iT}(2:end);
    D.b2.xsdt{iT} = xsdt2{iT};
    csdt2{iT} = csdt2{iT}(2:end);
    D.b2.csdt{iT} = csdt2{iT};    
    ih=1;
    for ir = ind_r
        % augment single replicats to obtain a cdf value for all x 
        csdr2{iT}(:,ih) = augment_cdf(xsd2{ir},xsdt2{iT},csd2{ir});
        % compute area between pooled and single replicat cdf
        D.areat2{iT}(ih) = sum(diff(xsdt2{iT}).*...
            abs(csdt2{iT}(1:end-1)-csdr2{iT}(1:end-1,ih)));
        ih=ih+1;
    end
    D.dist.b2.mean(iT) = mean(D.areat2{iT});
    D.dist.b2.var(iT) = var(D.areat2{iT})/length(ind_r);
end

% save data
save(options.name,'D','sol','options','parameters_true')
end