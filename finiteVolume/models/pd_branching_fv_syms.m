function [model] = pd_branching_fv_syms( )
% finite volume implementation for the pseudodynamics model with branching for Maehr data

model.param = 'log';
model.forward = true;
model.adjoint = false;

%% PDE discretization
% number of grid points on branch 1
n_grid1 = 300;

% number of grid points on branch 2
n_grid2 = 281;

grid = linspace(0,1,n_grid1);

% for length of grid interval, h, compute (1/h)^2
h2inv = (1/(grid(2)-grid(1)))^2;

% grid of centers of the grid intervals
grid_x = linspace(grid(1)+(grid(2)-grid(1))/2,grid(end-1)+(grid(end)-grid(end-1))/2,n_grid1-1); 

% branching region
grid_r = 20:123;

%% STATES
% create state syms: vector for x at each grid point
x = sym('x',[n_grid1+n_grid2-2,1]);

%% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms (values at spline nodes, 1-9 on branch 1; 10-12 on branch 2)
syms D1 D2 D3 D4 D5 D6 D7 D8 D9 D10 D11 D12 v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 delta1_2 delta2_1

p = [D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,delta1_2,delta2_1];

% compute symbolic value at grid

[Db1,~] = simpleSplineSym((0:0.125:1)',log([D1;D2;D3;D4;D5;D6;D7;D8;D9]),grid(2:end-1));
Db1 = exp(Db1);

% splines on grid(20):grid(end) => Interval 0:280/300
[Db2,~] = simpleSplineSym((0:0.5*(n_grid2-1)/(n_grid1-1):(n_grid2-1)/(n_grid1-1))',log([D10;D11;D12]),grid(2:n_grid2-1));
Db2 = exp(Db2);

[vb1,dvb1dgrid] = simpleSplineSym((0:0.125:1)',log([v1;v2;v3;v4;v5;v6;v7;v8;v9]),grid(2:end-1));
vb1 = exp(vb1);
dvb1dgrid = vb1.*dvb1dgrid;

% reduce v to 0 in the interval grid(end-29:end)
grid0 = 29;
[v_end_b1,~] = cspline_00(vb1(end-grid0),dvb1dgrid(end-grid0),grid(end-grid0),grid(end-1),grid(end-grid0:end-1));
vb1(end-grid0:end) = [v_end_b1,zeros(1,1)];

[vb2,dvb2dgrid] = simpleSplineSym((0:0.5*(n_grid2-1)/(n_grid1-1):(n_grid2-1)/(n_grid1-1))',log([v10;v11;v12]),grid(2:n_grid2-1));
vb2 = exp(vb2);
dvb2dgrid = vb2.*dvb2dgrid;

% reduce v to 0 in the interval grid(end-29:end)
grid0 = 29;
[v_end_b2,~] = cspline_00(vb2(end-grid0),dvb2dgrid(end-grid0),grid(n_grid2-grid0),grid(n_grid2-1),grid(n_grid2-grid0:n_grid2-1));
vb2(end-grid0:end) = [v_end_b2,zeros(1,1)];
                                  
ab1 = simpleSplineSym((0:0.125:1)',log([a1;a2;a3;a4;a5;a6;a7;a8;a9]),grid_x);
ab2 = simpleSplineSym((0:0.5*(n_grid2-1)/(n_grid1-1):(n_grid2-1)/(n_grid1-1))',log([a10;a11;a12]),grid_x(1:n_grid2-1));

%% CONSTANTS ( for these no sensitivities will be computed )
% this part is optional and can be ommited

% create parameter syms
k = sym('k',[n_grid1+n_grid2-2,1]);

%% SYSTEM EQUATIONS

% create symbolic variable for time
syms t

xdot = sym(zeros(size(x)));
% xdot is finite volume approximation to the PDE on the volumes defined by grid (n_grid-1 state variables)
% Robin boundary on the left hand side and Dirichlet on the right hand side

% dynamic for main branch
xdot(1) = h2inv*(-Db1(1)*(x(1)-x(2))) - vb1(1)*1/2*(x(1)+x(2))*sqrt(h2inv) + ab1(1)*x(1);

xdot(n_grid1-1) = h2inv*(Db1(end)*(x(n_grid1-2)-x(n_grid1-1)))+vb1(end)*1/2*(x(n_grid1-2)+x(n_grid1-1))*sqrt(h2inv) + ab1(n_grid1-1)*x(n_grid1-1);% 

for i = 2:grid_r(1)-1
 xdot(i) = h2inv*(Db1(i-1)*(x(i-1)-x(i))-Db1(i)*(x(i)-x(i+1))) + sqrt(h2inv)*1/2*(vb1(i-1)*(x(i-1)+x(i))-vb1(i)*(x(i)+x(i+1))) + ab1(i)*x(i);
end

% grid region
for i = grid_r
 xdot(i) = h2inv*(Db1(i-1)*(x(i-1)-x(i))-Db1(i)*(x(i)-x(i+1))) + sqrt(h2inv)*1/2*(vb1(i-1)*(x(i-1)+x(i))-vb1(i)*(x(i)+x(i+1))) + ab1(i)*x(i) - delta1_2*x(i) + delta2_1*x(n_grid1 + i - grid_r(1));
end

for i = grid_r(end)+1:n_grid1-2
 xdot(i) = h2inv*(Db1(i-1)*(x(i-1)-x(i))-Db1(i)*(x(i)-x(i+1))) + sqrt(h2inv)*1/2*(vb1(i-1)*(x(i-1)+x(i))-vb1(i)*(x(i)+x(i+1))) + ab1(i)*x(i);
end

% dynamic for side branch (from n_grid:1.8*n_grid-2)
xdot(n_grid1) = h2inv*(-Db2(1)*(x(n_grid1)-x(n_grid1+1))) - vb2(1)*1/2*(x(n_grid1)+x(n_grid1+1))*sqrt(h2inv) + ab2(1)*x(n_grid1) + delta1_2*x(grid_r(1)) - delta2_1*x(n_grid1);

xdot(n_grid1+n_grid2-2) = h2inv*(Db2(end)*(x(n_grid1+n_grid2-3)-x(n_grid1+n_grid2-2)))+vb2(end)*1/2*(x(n_grid1+n_grid2-3)+x(n_grid1+n_grid2-2))*sqrt(h2inv) + ab2(n_grid2-1)*x(n_grid1+n_grid2-2);% 

% grid region
for i = n_grid1+1:n_grid1+length(grid_r)-1
 xdot(i) = h2inv*(Db2(i-n_grid1)*(x(i-1)-x(i))-Db2(i-n_grid1+1)*(x(i)-x(i+1))) + sqrt(h2inv)*1/2*(vb2(i-n_grid1)*(x(i-1)+x(i))-vb2(i-n_grid1+1)*(x(i)+x(i+1))) + ab2(i-n_grid1+1)*x(i) + delta1_2*x(grid_r(i-n_grid1+1)) - delta2_1*x(i);
end

for i = n_grid1+length(grid_r):n_grid1+n_grid2-3
 xdot(i) = h2inv*(Db2(i-n_grid1)*(x(i-1)-x(i))-Db2(i-n_grid1+1)*(x(i)-x(i+1))) + sqrt(h2inv)*1/2*(vb2(i-n_grid1)*(x(i-1)+x(i))-vb2(i-n_grid1+1)*(x(i)+x(i+1))) + ab2(i-n_grid1+1)*x(i);
end

%% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

for i = 1:n_grid1+n_grid2-2
    x0(i) = k(i);
end

%% OBSERVALES
% number of cells on branch 1
N1 = sum(x(1:n_grid1-1)/sqrt(h2inv));
% number of cells on branch 2
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

