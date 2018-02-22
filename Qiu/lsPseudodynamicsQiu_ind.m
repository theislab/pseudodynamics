function [lS] = lsPseudodynamicsQiu_ind(varargin)
% compute the individual contribution to the  loglikelihood for Qiu data
% for the different data types

if nargin >= 3
	theta = varargin{1};
	modelFun = varargin{2};
	D = varargin{3};
else
	error('Not enough inputs!');
end

% default options
n_grid = 300;
options.grid = linspace(0,1,n_grid);
options.alpha = 100;
options.n_exp = 1;
if nargin == 4
	options = setdefault(varargin{4},options);
end

% u0 approximated by kernel density estimate
if ~isfield(options,'u0')
    options.u0 = ksdensity(D.ind.hist{1},x,'support',[0,1],'function','pdf');
    options.u0 = options.u0/trapz(x,options.u0)*exp(D.pop.mean(1));
    % compute initial cells in finite volumes
    options.u0 = 0.5*(options.u0(1:n_grid-1)+options.u0(2:n_grid));
else
    u0 = options.u0;
end

options_sim.sensi = 1;
options_sim.maxsteps = 1e6;
options_sim.sx0 = zeros(n_grid-1,length(theta));
% observation times
t = unique([D.pop.t;D.ind.t]);
% indices of pop size measurements in t
for iHelp = 1:length(D.pop.t)
    indP(iHelp) = find(t==D.pop.t(iHelp));
end
% indices of single cell measurements in t
for iHelp = 1:length(D.ind.t)
    indI(iHelp) = find(t==D.ind.t(iHelp));
end

% forward simulation
[status,t_sim,~,u,~,us] = modelFun(t,theta,u0,[],options_sim);

% variable asignment
N = u(:,end);
u = u(:,1:n_grid-1);

% computation of probability
p = bsxfun(@rdivide,u,N)*(1/(n_grid-1));

% computation of cdf
cs = cumsum(p');

% least-squares objective function
lS = zeros(1,3);

k=1;
for it = indI(2:end)
    x_combined = options.x_combined{k+1};
    Aug_matrix = options.Aug_matrix{k+1};
    cs_a = Aug_matrix*cs(:,it);
    csd_a = D.csd_a{k+1};
    I(k) = sum(diff(x_combined).*abs(cs_a(1:end-1)-csd_a(1:end-1)))^2;
    k=k+1;
end

% calculate MLE for sigma^2 
sigma2 = sum(I)/length(indI(2:end));

% evaluate negative log-likelihood for Least squares ansatz
for it = 1:length(indI)-1
    lS(1) = lS(1) + 0.5*log(2*pi*sigma2) + 1;
end
% log-likelihood value
lS(1) = -lS(1);

% negative log-likelihood for log pop size
for it = indP(2:end)
    lS(2) = lS(2) + 0.5*log(2*pi*D.pop.var(D.pop.t==t(it))/options.n_exp) + ...
        0.5*(D.pop.mean(D.pop.t==t(it))-log(N(it))).^2/...
        (D.pop.var(D.pop.t==t(it))/options.n_exp);
end
% log-likelihood value
lS(2) = -lS(2);

% regularization
lS(3) = lS(3) + options.alpha*(sum((theta(2:9)-theta(1:8)).^2)+...
    sum((theta(11:18)-theta(10:17)).^2)+sum((theta(20:27)-theta(19:26)).^2));
lS(3) = -lS(3);
end