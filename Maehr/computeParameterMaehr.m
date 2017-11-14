function [par] = computeParameterMaehr(parameters)
D1 = parameters(1:9);
D2 = parameters(10:12);
v1 = parameters(13:21);
v2 = parameters(22:24);
a1 = parameters(25:33);
a2 = parameters(34:36);

n_grid1 = 300;
n_grid2 = 281;
grid = linspace(0,1,n_grid1);

grid_x = linspace(grid(1)+(grid(2)-grid(1))/2,grid(end-1)+(grid(end)-grid(end-1))/2,n_grid1-1); 

grid_r = 20:123;
[Db1,~] = simpleSplineSym((0:0.125:1)',D1(:),grid(2:end-1));
par.Db1 = double(exp(Db1));
par.gridb1 = grid(2:end-1);
par.nodesb1 = 0:0.125:1;

[Db2,~] = simpleSplineSym((0:0.5*(n_grid2-1)/(n_grid1-1):(n_grid2-1)/(n_grid1-1))',D2(:),grid(2:n_grid2-1));
par.Db2 = double(exp(Db2));
par.gridb2 = grid(grid_r(1)+1:end-1);
par.nodesb2 = grid_r(1)+(0:0.5*(n_grid2-1)/(n_grid1-1):(n_grid2-1)/(n_grid1-1))';

grid0 = 29;
[vb1,dvb1dgrid] = simpleSplineSym((0:0.125:1)',v1(:),grid(2:end-1));
vb1 = exp(vb1);
dvb1dgrid = vb1.*dvb1dgrid;

grid0 = 29;
[v_end_b1,~] = cspline_00(vb1(end-grid0),dvb1dgrid(end-grid0),grid(end-grid0),grid(end-1),grid(end-grid0:end-1));
vb1(end-grid0:end) = [v_end_b1,zeros(1,1)];
par.vb1 = double(vb1);

[vb2,dvb2dgrid] = simpleSplineSym((0:0.5*(n_grid2-1)/(n_grid1-1):(n_grid2-1)/(n_grid1-1))',v2(:),grid(2:n_grid2-1));
vb2 = exp(vb2);
dvb2dgrid = vb2.*dvb2dgrid;

grid0 = 29;
[v_end_b2,~] = cspline_00(vb2(end-grid0),dvb2dgrid(end-grid0),grid(n_grid2-grid0),grid(n_grid2-1),grid(n_grid2-grid0:n_grid2-1));
vb2(end-grid0:end) = [v_end_b2,zeros(1,1)];
par.vb2 = double(vb2);

par.ab1 = double(simpleSplineSym((0:0.125:1)',a1(:),grid_x));
par.ab2 = double(simpleSplineSym((0:0.5*(n_grid2-1)/(n_grid1-1):(n_grid2-1)/(n_grid1-1))',a2(:),grid_x(1:n_grid2-1)));

par.gridxb1 = grid_x;
par.gridxb2 = grid_x(grid_r(1):end);

par.d12 = exp(parameters(end-1));
par.d21 = exp(parameters(end));
end