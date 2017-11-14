function [lS,grad] = lsPseudodynamicsKlein(varargin)
% compute the least-squares difference between measured and simulated area

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
    u0 = ksdensity(D.ind.hist{1},options.grid,'support',[0,1],'function','pdf');
    u0 = 0.5*([u0(1:n_grid-1)+u0(2:n_grid),u0(n_grid+1:end-1)+u0(n_grid+2:end)]);
else
    u0 = options.u0;
end

options_sim.sensi = 1;
options_sim.maxsteps = 1e6;
options_sim.sx0 = zeros(n_grid-1,length(theta));
% forward simulation
[status,t_sim,~,u,~,us] = modelFun(D.pop.t,theta,u0,[],options_sim);

% variable asignment
N = u(:,end);
dNdtheta = us(:,end,:);

u = u(:,1:n_grid-1);
dudtheta = us(:,1:n_grid-1,:);

% computation of probability
p = bsxfun(@rdivide,u,N)*(1/(n_grid-1));
dpdtheta = 1/(n_grid-1)*bsxfun(@rdivide,bsxfun(@times,N,dudtheta)-bsxfun(@times,u,dNdtheta),N.^2);

% computation of cdf
cs = cumsum(p');
dcsdtheta = cumsum(dpdtheta,2);

% least-squares objective function
lS = 0;
grad = zeros(size(theta));
        
area = nan(length(D.pop.t)-1,1);
dareadtheta = nan(length(D.pop.t)-1,length(theta));
     
for it = 2:length(D.pop.t)
    x_combined = options.x_combined{it};
    Aug_matrix = options.Aug_matrix{it};
    % augment simulated cdf
    cs_a = Aug_matrix*cs(:,it);
    dcs_adtheta = Aug_matrix*squeeze(dcsdtheta(it,:,:));
    csd_a = D.csd_a{it};
    lS = lS + sum(diff(x_combined).*abs(cs_a(1:end-1)-csd_a(1:end-1)))^2;
    grad = grad + (2*sum(diff(x_combined).*abs(cs_a(1:end-1)-csd_a(1:end-1)))*(diff(x_combined).*sign(cs_a(1:end-1)-csd_a(1:end-1)))'*(dcs_adtheta(1:end-1,:)))';
end

% regularization
lS = lS + options.alpha*(sum((theta(2:9)-theta(1:8)).^2)+sum((theta(11:18)-theta(10:17)).^2));
grad = grad + options.alpha*2*([-(theta(2:9)-theta(1:8));0;-(theta(11:18)-theta(10:17));0]...
    +[0;theta(2:9)-theta(1:8);0;theta(11:18)-theta(10:17)]);
end