function [varargout] = computeSolForward(parameters,D,options)
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

%grid_r = 20:123;
[Db1,~] = simpleSplineSym((0:0.125:1)',D1(:),grid(2:end-1));
Db1 = exp(Db1);

[Db2,~] = simpleSplineSym((0:0.5*(n_grid2-1)/(n_grid1-1):(n_grid2-1)/(n_grid1-1))',D2(:),grid(2:n_grid2-1));
Db2 = exp(Db2);

[vb1,dvb1dgrid] = simpleSplineSym((0:0.125:1)',v1(:),grid(2:end-1));
vb1 = exp(vb1);
dvb1dgrid = vb1.*dvb1dgrid;

grid0 = 29;
[v_end_b1,~] = cspline_00(vb1(end-grid0),dvb1dgrid(end-grid0),grid(end-grid0),grid(end-1),grid(end-grid0:end-1));
vb1(end-grid0:end) = [v_end_b1,zeros(1,1)];

[vb2,dvb2dgrid] = simpleSplineSym((0:0.5*(n_grid2-1)/(n_grid1-1):(n_grid2-1)/(n_grid1-1))',v2(:),grid(2:n_grid2-1));
vb2 = exp(vb2);
dvb2dgrid = vb2.*dvb2dgrid;

grid0 = 29;
[v_end_b2,~] = cspline_00(vb2(end-grid0),dvb2dgrid(end-grid0),grid(n_grid2-grid0),grid(n_grid2-1),grid(n_grid2-grid0:n_grid2-1));
vb2(end-grid0:end) = [v_end_b2,zeros(1,1)];

ab1 = simpleSplineSym((0:0.125:1)',a1(:),grid_x);
ab2 = simpleSplineSym((0:0.5*(n_grid2-1)/(n_grid1-1):(n_grid2-1)/(n_grid1-1))',a2(:),grid_x(1:n_grid2-1));

theta = [Db1,Db2,vb1,vb2,ab1,ab2,exp(parameters(end-1)),exp(parameters(end))];
theta = double(theta);

sol = simulate_pd_branching_forward(D.pop.t,theta,options.u0,[],[]);
sol.p1 = bsxfun(@rdivide,sol.y(:,1:299),sol.y(:,end-1));
sol.p2 = bsxfun(@rdivide,sol.y(:,300:end-2),sol.y(:,end));

if nargout == 1
    varargout{1} = sol;
elseif nargout == 2;
    varargout{1} = sol;
    varargout{2} = theta;
end
end

