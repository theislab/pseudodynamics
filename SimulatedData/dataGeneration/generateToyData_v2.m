function [] = generateToyData_v2(varargin)
% set model specifications
n_grid = 300;
x = linspace(0,1,n_grid);

grid_x = linspace(x(1)+(x(2)-x(1))/2,x(end-1)+(x(end)-x(end-1))/2,n_grid-1); 

% assign inputs
parameters_true = varargin{1};
modelfun = varargin{2};
D = varargin{3};
options = varargin{4};

% simulate model
sol = modelfun(D.pop.t,parameters_true,options.u0);
sol.p = bsxfun(@rdivide,sol.y(:,1:end-1),sol.y(:,end))*(1/(n_grid-1));

% generate count data from distribution
cdf_sim = cumsum(sol.p,2);
k=1;
for iTime = D.ind.tp
    [cdf,mask] = unique(cdf_sim(D.pop.t == iTime,:));
    xm = x(mask);
    % generate uniform random number
    rS = rand(1,options.nSample);
    % draw uniform number in interval
    for iS=1:options.nSample
        lb = xm(find(cdf>rS(iS),1,'first'));
        D.ind.hist{k}(iS) = lb+1/299*rand(1);
    end
    D.ind.cdf{k} = histcounts(D.ind.hist{k},x,'Normalization','cdf');
    k=k+1;
end

for iT=1:length(D.pop.t)
    D.ind.size{iT} = normrnd(sol.y(iT,end),sqrt(options.varPop(iT)),[1,options.pSample]);
end
% compute area statistics
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
    %D.pop.mean(it) = mean(D.ind.size(ind_r));
    %D.pop.var(it) = var(D.ind.size(ind_r))/length(ind_r);
    D.pop.mean(it) = mean(D.ind.size{it});
    D.pop.var(it) = var(D.ind.size{it})/options.pSample;
    % pool all data for a time point
    for ir = ind_r
        D_t = [D_t,D.ind.hist{ir}];
    end
    % empirical cdf for pooled data
    D.dist.cdf{it} = histcounts(D_t,x,'Normalization','cdf');
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
        D.areat{it}(ih) = sum(diff(xsdt{it}).*abs(csdt{it}(1:end-1)-csdr{it}(1:end-1,ih)));
        D.areat_cdf{it}(ih) = sum(diff(x).*abs(D.ind.cdf{ir}-D.dist.cdf{it}));
        ih=ih+1;
    end
    x_aug = union(x(1:end-1),xsdt{it});
    D.csd_a{it} = augment_cdf(xsdt{it},x_aug,csdt{it});
    cdf_aug = augment_cdf(x(1:end-1),x_aug,cdf_sim(it,:)');
    D.area_true(it) = sum(diff(x_aug).*abs(D.csd_a{it}(1:end-1)-cdf_aug(1:end-1)));
    D.area_cdf_true(it) = sum(diff(x).*abs(D.dist.cdf{it}-cdf_sim(it,:)));
    D.dist.mean(it) = mean(D.areat{it});
    D.dist.mean_cdf(it) = mean(D.areat_cdf{it});
    D.dist.var(it) = var(D.areat{it})/length(ind_r);
    D.dist.var_cdf(it) = var(D.areat_cdf{it})/length(ind_r);
end

figure; plot(D.area_true)
hold on; plot(D.area_cdf_true)
hold on; plot(D.dist.mean,'bo')
hold on; plot(D.dist.mean_cdf,'ro')
hold on; plot(D.dist.mean+2*sqrt(D.dist.var),'bx')
hold on; plot(D.dist.mean-2*sqrt(D.dist.var),'bx')
hold on; plot(D.dist.mean_cdf+2*sqrt(D.dist.var_cdf),'rx')
hold on; plot(D.dist.mean_cdf-2*sqrt(D.dist.var_cdf),'rx')
legend('true - augmented grid','true - fixed grid','mean - ag','mean - fg','2 sigma - ag','','2 sigma - fg','')
% save data
save(options.name,'D','sol','options','parameters_true')
end
