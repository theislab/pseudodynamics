function varargout = lsPseudodynamicsFA(varargin)
% compute the least-squares difference to the KO data, only on branch 1

if nargin >= 2
	modelFun = varargin{1};
	D = varargin{2};
else
	error('Not enough inputs!');
end

n_grid = 300;
% default options
options.grid = linspace(0,1,n_grid);
options.alpha = 100;
options.n_exp = 1;
if nargin == 3
	options = setdefault(varargin{3},options);
end

% forward simulation
sol = modelFun;
u = sol.y;

N1 = u(:,end-1);
u1 = u(:,1:n_grid-1);

% computation of probability
p1 = bsxfun(@rdivide,u1,N1)*(1/(n_grid-1));
p_plot = bsxfun(@rdivide,u1,N1); 

% figure; 
% hold on;
% for i=1:length(D.ind.Rag.hist)
% histogram(D.ind.Rag.hist{i},'normalization','pdf','DisplayStyle','stairs','Linewidth',2);
% end
% stairs(options.grid(1:end-1),p_plot(:,1:299)','Linewidth',2)
% xlabel('dpt')
% ylabel('probability')
% legend('data1','data2','data3','sim')

% compute simulated cdf
cs1 = cumsum(p1');

lS = 0;
for it = 2:length(D.pop.t)
    x_combined = options.x_combined{it};
    Aug_matrix = options.Aug_matrix{it};
    cs_a = Aug_matrix*cs1(:,it);
    csd_a = D.csd_a{it};
    lS = lS + sum(diff(x_combined).*abs(cs_a(1:end-1)-csd_a(1:end-1)))^2;
end

if nargout==1
    varargout{1}=lS;
elseif nargout == 2;
    varargout{1} = lS;
    varargout{2} = p_plot;
end