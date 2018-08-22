function [negLogL,grad] = llCompartmentQiut4(varargin)
% compute the log likelihood value for Qiu compartment data for
% time-dependent rates

if nargin >= 3
	theta = varargin{1};
	modelFun = varargin{2};
	D = varargin{3};
else
	error('Not enough inputs!');
end

% default options
options.n_exp = 1;
options.alpha = 0;
if nargin == 4
	options = setdefault(varargin{4},options);
end

u0 = options.u0;

options_sim.sensi = 1;
options_sim.maxsteps = 1e6;
options_sim.sx0 = zeros(2,length(theta));
% observation times
t = unique([D.pop.t;D.ind.tp]);
% indices of pop size measurements in t
for iHelp = 1:length(D.pop.t)
    indP(iHelp) = find(t==D.pop.t(iHelp));
end
% indices of single cell measurements in t
for iHelp = 1:length(D.ind.tp(2:end))
    indI(iHelp) = find(t==D.ind.t(iHelp));
end

% forward simulation
[status,~,~,u,~,us] = modelFun(t,theta,u0,[],options_sim);
% check if simulation failed
if status ~= 0
    error([],'AMICI simulation failed')
end
% check if solution is negative
if max(max(u<0))
    error('negative solution')
end

% variable asignment
N = u(:,end);
dNdtheta = us(:,end,:);

% mature in 2nd compartment
u = u(:,2);
dudtheta = us(:,2,:);

% computation of probability
p = bsxfun(@rdivide,u,N);
dpdtheta = bsxfun(@rdivide,bsxfun(@times,N,dudtheta)-bsxfun(@times,u,dNdtheta),N.^2);

% log-likelihood
negLogL = 0;
grad = zeros(size(theta));

% negative log-likelihood for fraction in mature compartment
for it = indI
     negLogL = negLogL + 0.5*log(2*pi*D.mature.var(D.ind.t==t(it))) +...
         0.5*(D.mature.mean(D.ind.t==t(it))-p(it)).^2/(D.mature.var(D.ind.t==t(it)));
    grad = grad  - ((D.mature.mean(D.ind.t==t(it))-p(it))/...
        (D.mature.var(D.ind.t==t(it)))*squeeze(dpdtheta(it,:))');
end

% negative log-likelihood for log pop size
for it = indP(2:end)    
    negLogL = negLogL + 0.5*log(2*pi*D.pop.var(D.pop.t==t(it))/options.n_exp) + ...
        0.5*(D.pop.mean(D.pop.t==t(it))-log(N(it))).^2/...
        (D.pop.var(D.pop.t==t(it))/options.n_exp);
    grad = grad  - ((D.pop.mean(D.pop.t == t(it))-log(N(it)))/...
        (D.pop.var(D.pop.t == t(it))/options.n_exp)*1/(N(it))*squeeze(dNdtheta(it,:,:))')';
end

% regularization
negLogL = negLogL + options.alpha*(sum((theta(2:3)-theta(1:2)).^2)+...
    sum((theta(5:6)-theta(4:5)).^2)+sum((theta(8:9)-theta(7:8)).^2)+...
    sum((theta(11:12)-theta(10:11)).^2));
grad = grad + options.alpha*2*([-(theta(2:3)-theta(1:2));0;...
    -(theta(5:6)-theta(4:5));0;-(theta(8:9)-theta(7:8));0;...
    -(theta(11:12)-theta(10:11));0]...
    +[0;theta(2:3)-theta(1:2);0;theta(5:6)-theta(4:5);0;...
    theta(8:9)-theta(7:8);0;theta(11:12)-theta(10:11)]);
end