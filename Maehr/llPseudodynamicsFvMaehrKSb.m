function [LogL,grad] = llPseudodynamicsFvMaehrKSb(varargin)
% compute likelihood for pseudodynamics

if nargin >= 3
	theta = varargin{1};
	modelFun = varargin{2};
	D = varargin{3};
else
	error('Not enough inputs!');
end

n_grid = 300;
grid_r = 20:123;
% default options
options.grid = linspace(0,1,n_grid);
options.alpha = 100;
options.n_exp = 1;
if nargin == 4
	options = setdefault(varargin{4},options);
end

% u0 approximated by kernel density estimate
if ~isfield(options,'u0')
    % compute mean u0 distribution
    for i0 = 1:3
        u0(i0,:) = ksdensity(D.ind.b1.hist{i0},x,'support',[0,1],'function','pdf');
        N01(i0) = length(D.ind.b1.hist{i0});
        u0(i0,:) = u0(i0,:)/trapz(x,u0(i0,:));%*D.ind.size(i0);
    end
    for i0 = 1:3
        u02(i0,:) = ksdensity(D.ind.b2.hist{i0},x(20:end),'support',[0,1],'function','pdf');
        N02(i0) = length(D.ind.b2.hist{i0});
        u02(i0,:) = u02(i0,:)/trapz(x(20:end),u02(i0,:));%*D.ind.size(i0);
    end
    % weighting of subpopulations
    w1=mean(N01./(N01+N02));
    w2=mean(N02./(N01+N02));

    u0 = [w1*mean(u0,1),w2*mean(u02,1)]*D.pop.mean(1);

    u0 = 0.5*([u0(1:n_grid-1)+u0(2:n_grid),u0(n_grid+1:end-1)+u0(n_grid+2:end)]);
else
    u0 = options.u0;
end

options_sim.sensi = 1;
options_sim.maxsteps = 1e6;
options_sim.sx0 = zeros(n_grid+n_grid-19-2,length(theta));
% simulation
[status,t_sim,~,u,~,us] = modelFun(D.pop.t,theta,u0,[],options_sim);

% variable assignment
N1 = u(:,end-1);
dN1dtheta = us(:,end-1,:);

N2 = u(:,end);
dN2dtheta = us(:,end,:);

N = N1+N2;
dNdtheta = dN1dtheta + dN2dtheta;

u1 = u(:,1:n_grid-1);
du1dtheta = us(:,1:n_grid-1,:);

u2 = u(:,n_grid:end-2);
du2dtheta = us(:,n_grid:end-2,:);

% computation of probability on individual branches
p1 = bsxfun(@rdivide,u1,N1)*(1/(n_grid-1));
dp1dtheta = 1/(n_grid-1)*bsxfun(@rdivide,bsxfun(@times,N1,du1dtheta)-bsxfun(@times,u1,dN1dtheta),N1.^2);

p2 = bsxfun(@rdivide,u2,N2)*(1/(n_grid-1));
dp2dtheta = 1/(n_grid-1)*bsxfun(@rdivide,bsxfun(@times,N2,du2dtheta)-bsxfun(@times,u2,dN2dtheta),N2.^2);

% compute simulated cdf
cs1 = cumsum(p1');
dcs1dtheta = cumsum(dp1dtheta,2);

cs2 = cumsum(p2');
dcs2dtheta = cumsum(dp2dtheta,2);

% compute simulated weights
w1 = N1./N;
dw1dtheta = bsxfun(@rdivide,bsxfun(@times,N,dN1dtheta)-bsxfun(@times,N1,dNdtheta),N.^2);

% logLikelihood
LogL = 0;
grad = zeros(size(theta));

area1 = nan(length(D.pop.t)-1,1);
darea1dtheta = nan(length(D.pop.t)-1,length(theta));
area2 = nan(length(D.pop.t),1);
darea2dtheta = nan(length(D.pop.t),length(theta));

% area statistics branch 1
for it = 2:length(D.pop.t)
    x_combined = options.x_combined{1,it};
    Aug_matrix = options.Aug_matrix{1,it};
    cs1_a = Aug_matrix*cs1(:,it);
    dcs1_adtheta = Aug_matrix*squeeze(dcs1dtheta(it,:,:));
    csd1_a = D.b1.csd_a{it};
    area1(it) = sum(diff(x_combined).*abs(cs1_a(1:end-1)-csd1_a(1:end-1)));
    darea1dtheta(it,:) = (diff(x_combined).*sign(cs1_a(1:end-1)-csd1_a(1:end-1)))'*(dcs1_adtheta(1:end-1,:));
    LogL = LogL + 0.5*log(2*pi*D.dist.b1.var(it)) + 0.5*(D.dist.b1.mean(it)-area1(it)).^2/D.dist.b1.var(it);
    grad = grad - ((D.dist.b1.mean(it)-area1(it))/D.dist.b1.var(it)*darea1dtheta(it,:))';
end

% area statistics branch 2
for it = 2:length(D.pop.t)
    x_combined = options.x_combined{2,it};
    Aug_matrix = options.Aug_matrix{2,it};
    cs2_a = Aug_matrix*cs2(:,it);
    dcs2_adtheta = Aug_matrix*squeeze(dcs2dtheta(it,:,:));
    csd2_a = D.b2.csd_a{it};
    area2(it) = sum(diff(x_combined).*abs(cs2_a(1:end-1)-csd2_a(1:end-1)));
    darea2dtheta(it,:) = (diff(x_combined).*sign(cs2_a(1:end-1)-csd2_a(1:end-1)))'*(dcs2_adtheta(1:end-1,:));
    LogL = LogL +  0.5*log(2*pi*D.dist.b2.var(it)) + 0.5*(D.dist.b2.mean(it)-area2(it)).^2/D.dist.b2.var(it);
    grad = grad  - ((D.dist.b2.mean(it)-area2(it))/D.dist.b2.var(it)*darea2dtheta(it,:))';
end

% population level statistics
for it = 2:length(D.pop.t)
    LogL = LogL + 0.5*log(2*pi*D.pop.var(it)/options.n_exp) + 0.5*(D.pop.mean(it)-N(it)).^2/(D.pop.var(it)/options.n_exp);
    grad = grad  - ((D.pop.mean(it)-N(it))/(D.pop.var(it)/options.n_exp)*squeeze(dNdtheta(it,:,:))')';
                    
    LogL = LogL + 0.5*log(2*pi*D.weight.b1.var(it)) + 0.5*(D.weight.b1.mean(it)-w1(it)).^2/(D.weight.b1.var(it));
    grad = grad  - ((D.weight.b1.mean(it)-w1(it))/(D.weight.b1.var(it))*squeeze(dw1dtheta(it,:,:))')';
end

% regularization
LogL = LogL + options.alpha*(sum((theta(2:9)-theta(1:8)).^2)+sum((theta(11:12)-theta(10:11)).^2)+sum((theta(14:21)-theta(13:20)).^2)+sum((theta(23:24)-theta(22:23)).^2)+sum((theta(26:33)-theta(25:32)).^2)+sum((theta(35:36)-theta(34:35)).^2));
grad = grad + options.alpha*2*([-(theta(2:9)-theta(1:8));0;-(theta(11:12)-theta(10:11));0;-(theta(14:21)-theta(13:20));0;-(theta(23:24)-theta(22:23));0;-(theta(26:33)-theta(25:32));0;-(theta(35:36)-theta(34:35));0;0;0]...
    +[0;theta(2:9)-theta(1:8);0;theta(11:12)-theta(10:11);0;theta(14:21)-theta(13:20);0;theta(23:24)-theta(22:23);0;theta(26:33)-theta(25:32);0;theta(35:36)-theta(34:35);0;0]);
end