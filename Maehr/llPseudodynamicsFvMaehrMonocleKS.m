function [LogL,grad] = llPseudodynamicsFvMaehrMonocleKS(varargin)
% compute likelihood for pseudodynamics on Monocle processed Maehr data

if nargin >= 3
	theta = varargin{1};
	modelFun = varargin{2};
	D = varargin{3};
else
	error('Not enough inputs!');
end

n_grid = 300;
% default options
options.grid = linspace(0,1,n_grid);
options.alpha = 100;
options.n_exp = 1;
if nargin == 4
	options = setdefault(varargin{4},options);
end

% u0 approximated by kernel density estimate
if ~isfield(options,'u0')
    for i0 = 1:3
        u0(i0,:) = ksdensity(D.ind.hist{i0},x,'support',[0,1],'function','pdf');
        N01(i0) = length(D.ind.hist{i0});
        u0(i0,:) = u0(i0,:)/trapz(x,u0(i0,:));
    end
    u0 = mean(u0,1)*D.pop.mean(1);
    u0 = 0.5*(u0(1:n_grid-1)+u0(2:n_grid));
else
    u0 = options.u0;
end

options_sim.sensi = 1;
options_sim.maxsteps = 1e6;
options_sim.sx0 = zeros(n_grid-1,length(theta));
% simulation
[status,t_sim,~,u,~,us] = modelFun(D.pop.t,theta,u0,[],options_sim);

% variable assignment
N = u(:,end);
dNdtheta = us(:,end,:);

u = u(:,1:n_grid-1);
dudtheta = us(:,1:n_grid-1,:);


% computation of probability
p = bsxfun(@rdivide,u,N)*(1/(n_grid-1));
dpdtheta = 1/(n_grid-1)*bsxfun(@rdivide,bsxfun(@times,N,dudtheta)-bsxfun(@times,u,dNdtheta),N.^2);

% compute simulated cdf
cs = cumsum(p');
dcsdtheta = cumsum(dpdtheta,2);

% logLikelihood
LogL = 0;
grad = zeros(size(theta));

area = nan(length(D.pop.t)-1,1);
dareadtheta = nan(length(D.pop.t)-1,length(theta));

% area statistics
for it = 2:length(D.pop.t)
    x_combined = options.x_combined{it};
    Aug_matrix = options.Aug_matrix{it};
    cs_a = Aug_matrix*cs(:,it);
    dcs_adtheta = Aug_matrix*squeeze(dcsdtheta(it,:,:));
    csd_a = D.csd_a{it};
    area(it) = sum(diff(x_combined).*abs(cs_a(1:end-1)-csd_a(1:end-1)));
    dareadtheta(it,:) = (diff(x_combined).*sign(cs_a(1:end-1)-csd_a(1:end-1)))'*(dcs_adtheta(1:end-1,:));
    LogL = LogL + 0.5*log(2*pi*D.dist.var(it)) + 0.5*(D.dist.mean(it)-area(it)).^2/D.dist.var(it);
    grad = grad - ((D.dist.mean(it)-area(it))/D.dist.var(it)*dareadtheta(it,:))';
end

% population level statistics
for it = 2:length(D.pop.t)
    LogL = LogL + 0.5*log(2*pi*D.pop.var(it)/options.n_exp) + 0.5*(D.pop.mean(it)-N(it)).^2/(D.pop.var(it)/options.n_exp);
    grad = grad  - ((D.pop.mean(it)-N(it))/(D.pop.var(it)/options.n_exp)*squeeze(dNdtheta(it,:,:))')';
end

% regularization
LogL = LogL + options.alpha*(sum((theta(2:9)-theta(1:8)).^2)+sum((theta(11:18)-theta(10:17)).^2)+sum((theta(20:27)-theta(19:26)).^2));
grad = grad + options.alpha*2*([-(theta(2:9)-theta(1:8));0;-(theta(11:18)-theta(10:17));0;-(theta(20:27)-theta(19:26));0]...
    +[0;theta(2:9)-theta(1:8);0;theta(11:18)-theta(10:17);0;theta(20:27)-theta(19:26)]);
end