function [par] = computeParameterMaehrMonocle(parameters)
D = parameters(1:9);
v = parameters(10:18);
a = parameters(19:27);

n_grid1 = 300;
grid = linspace(0,1,n_grid1);

grid_x = linspace(grid(1)+(grid(2)-grid(1))/2,grid(end-1)+(grid(end)-grid(end-1))/2,n_grid1-1); 

[D,~] = simpleSplineSym((0:0.125:1)',D(:),grid(2:end-1));
par.D = double(exp(D));
par.grid = grid(2:end-1);
par.nodes = 0:0.125:1;

[v,dvdgrid] = simpleSplineSym((0:0.125:1)',v(:),grid(2:end-1));
v = exp(v);
dvdgrid = v.*dvdgrid;

grid0 = 29;
[v_end,~] = cspline_00(v(end-grid0),dvdgrid(end-grid0),grid(end-grid0),grid(end-1),grid(end-grid0:end-1));
v(end-grid0:end) = [v_end,zeros(1,1)];
par.v = double(v);


par.a = double(simpleSplineSym((0:0.125:1)',a(:),grid_x));

par.gridx = grid_x;

end