function [model] = pd_branching_Rag_fv_forward_syms( )
% finite volume implementation 

model.param = 'lin';
model.forward = false;
model.adjoint = false;

%% PDE discretization
n_grid1 = 300;
n_grid2 = 279;

grid = linspace(0,1,n_grid1);

% (1/h)^2
h2inv = (1/(grid(2)-grid(1)))^2;

grid_x = linspace(grid(1)+(grid(2)-grid(1))/2,grid(end-1)+(grid(end)-grid(end-1))/2,n_grid1-1); 

grid_r = 22:111;
%% STATES
% create state syms: vector for x at each grid point
x = sym('x',[n_grid1+n_grid2-2,1]);

%% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms
% number of volumes: n_grid1-1 + n_grid2-1
% number of non zero D and v entries: n_grid1-2 + n_grid2-2
D = sym('D',[1,n_grid1+n_grid2-4]);
v = sym('v',[1,n_grid1+n_grid2-4]);
a = sym('a',[1,n_grid1+n_grid2-2]);
syms delta1_2 delta2_1

%p = [D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,delta1_2,delta2_1];
p = [D,v,a,delta1_2,delta2_1];

%% CONSTANTS ( for these no sensitivities will be computed )
% this part is optional and can be ommited

% create parameter syms
k = sym('k',[n_grid1+n_grid2-2,1]);

%% SYSTEM EQUATIONS

% create symbolic variable for time
syms t

xdot = sym(zeros(size(x)));
% xdot is finite difference approximation to PDE on grid
% boundary x(0) = 0;

% dynamic for main branch
xdot(1) = h2inv*(-D(1)*(x(1)-x(2))) - v(1)*1/2*(x(1)+x(2))*sqrt(h2inv) + a(1)*x(1);

xdot(n_grid1-1) = h2inv*(D(n_grid1-2)*(x(n_grid1-2)-x(n_grid1-1)))+v(n_grid1-2)*1/2*(x(n_grid1-2)+x(n_grid1-1))*sqrt(h2inv) + a(n_grid1-1)*x(n_grid1-1);% 

for i = 2:grid_r(1)-1
 xdot(i) = h2inv*(D(i-1)*(x(i-1)-x(i))-D(i)*(x(i)-x(i+1))) + sqrt(h2inv)*1/2*(v(i-1)*(x(i-1)+x(i))-v(i)*(x(i)+x(i+1))) + a(i)*x(i);
end

for i = grid_r
 xdot(i) = h2inv*(D(i-1)*(x(i-1)-x(i))-D(i)*(x(i)-x(i+1))) + sqrt(h2inv)*1/2*(v(i-1)*(x(i-1)+x(i))-v(i)*(x(i)+x(i+1))) + a(i)*x(i) - delta1_2*x(i) + delta2_1*x(n_grid1 + i - grid_r(1));
end

for i = grid_r(end)+1:n_grid1-2
 xdot(i) = h2inv*(D(i-1)*(x(i-1)-x(i))-D(i)*(x(i)-x(i+1))) + sqrt(h2inv)*1/2*(v(i-1)*(x(i-1)+x(i))-v(i)*(x(i)+x(i+1))) + a(i)*x(i);
end

% dynamic for side branch (from n_grid:1.8*n_grid-2)
xdot(n_grid1) = h2inv*(-D(n_grid1-1)*(x(n_grid1)-x(n_grid1+1))) - v(n_grid1-1)*1/2*(x(n_grid1)+x(n_grid1+1))*sqrt(h2inv) + a(n_grid1)*x(n_grid1) + delta1_2*x(grid_r(1)) - delta2_1*x(n_grid1);

xdot(n_grid1+n_grid2-2) = h2inv*(D(end)*(x(n_grid1+n_grid2-3)-x(n_grid1+n_grid2-2)))+v(end)*1/2*(x(n_grid1+n_grid2-3)+x(n_grid1+n_grid2-2))*sqrt(h2inv) + a(n_grid1+n_grid2-2)*x(n_grid1+n_grid2-2);% 

for i = n_grid1+1:n_grid1+length(grid_r)-1
 xdot(i) = h2inv*(D(i-2)*(x(i-1)-x(i))-D(i-1)*(x(i)-x(i+1))) + sqrt(h2inv)*1/2*(v(i-2)*(x(i-1)+x(i))-v(i-1)*(x(i)+x(i+1))) + a(i)*x(i) + delta1_2*x(grid_r(i-n_grid1+1)) - delta2_1*x(i);
end

for i = n_grid1+length(grid_r):n_grid1+n_grid2-3
 xdot(i) = h2inv*(D(i-2)*(x(i-1)-x(i))-D(i-1)*(x(i)-x(i+1))) + sqrt(h2inv)*1/2*(v(i-2)*(x(i-1)+x(i))-v(i-1)*(x(i)+x(i+1))) + a(i)*x(i);
end

%% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

for i = 1:n_grid1+n_grid2-2
    x0(i) = k(i);
end

%% OBSERVALES
N1 = sum(x(1:n_grid1-1)/sqrt(h2inv));
N2 = sum(x(n_grid1:end)/sqrt(h2inv));
y = [x;N1;N2];


%% SYSTEM STRUCT

model.sym.x = x;
model.sym.k = k;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.x0 = x0;
model.sym.y = y;
end

