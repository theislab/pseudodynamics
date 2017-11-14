function [negLogL] = llPseudodynamicsFvMaehrCV(varargin)
% compute likelihood for corssvalidation prediction error for pseudodynamics
% input here is the coputed solution not the parameter

if nargin >= 3
	theta = varargin{1};
	sol = varargin{2};
	D = varargin{3};
else
	error('Not enough inputs!');
end

n_grid = 300;
% default options
options.alpha = 100;
options.n_exp = 1;
options.it = 2:length(D.pop.t);
if nargin == 4
	options = setdefault(varargin{4},options);
end

N1 = sol.y(:,end-1);
N2 = sol.y(:,end);

N = N1+N2;

u1 = sol.y(:,1:n_grid-1);
u2 = sol.y(:,n_grid:end-2);

% computation of probability
p1 = bsxfun(@rdivide,u1,N1)*(1/(n_grid-1));
p2 = bsxfun(@rdivide,u2,N2)*(1/(n_grid-1));

cs1 = cumsum(p1');
cs2 = cumsum(p2');

w1 = N1./N;

% logLikelihood without regularization only for time point it
negLogL = 0;

area1 = nan(length(options.it),1);
area2 = nan(length(D.pop.t)-1,1);

% area statistics branch 1
for it = options.it
    x_combined = options.x_combined{1,it};
    Aug_matrix = options.Aug_matrix{1,it};
    cs1_a = Aug_matrix*cs1(:,it);
    csd1_a = D.b1.csd_a{it};
    area1(it) = sum(diff(x_combined).*abs(cs1_a(1:end-1)-csd1_a(1:end-1)));
    negLogL = negLogL + 0.5*log(2*pi*D.dist.b1.var(it)) + 0.5*(D.dist.b1.mean(it)-area1(it)).^2/D.dist.b1.var(it);
end

% area statistics branch 2
for it = options.it
    x_combined = options.x_combined{2,it};
    Aug_matrix = options.Aug_matrix{2,it};
    cs2_a = Aug_matrix*cs2(:,it);    
    csd2_a = D.b2.csd_a{it};
    area2(it) = sum(diff(x_combined).*abs(cs2_a(1:end-1)-csd2_a(1:end-1)));
    negLogL = negLogL +  0.5*log(2*pi*D.dist.b2.var(it)) + 0.5*(D.dist.b2.mean(it)-area2(it)).^2/D.dist.b2.var(it);
end

% population level statistics
for it = options.it
    negLogL = negLogL + 0.5*log(2*pi*D.pop.var(it)/options.n_exp) + 0.5*(D.pop.mean(it)-N(it)).^2/(D.pop.var(it)/options.n_exp);
    negLogL = negLogL + 0.5*log(2*pi*D.weight.b1.var(it)) + 0.5*(D.weight.b1.mean(it)-w1(it)).^2/(D.weight.b1.var(it));
end

% regularization
% negLogL = negLogL + options.alpha*(sum((theta(2:9)-theta(1:8)).^2)+sum((theta(11:12)-theta(10:11)).^2)+sum((theta(14:21)-theta(13:20)).^2)+sum((theta(23:24)-theta(22:23)).^2)+sum((theta(26:33)-theta(25:32)).^2)+sum((theta(35:36)-theta(34:35)).^2));
end