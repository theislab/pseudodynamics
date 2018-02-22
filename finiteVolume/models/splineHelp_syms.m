function [model] = splineHelp_syms( )
% simulation file for the evaluation of the spline for the growth rate 
% at different timepoints 

model.param = 'log';
model.forward = true;
model.adjoint = false;

%% PDE discretization
n_grid = 300;
grid = linspace(0,1,n_grid);

h2inv = (1/(grid(2)-grid(1)))^2;

grid_x = linspace(grid(1)+(grid(2)-grid(1))/2,grid(end-1)+(grid(end)-grid(end-1))/2,n_grid-1); 
%% STATES
% create state syms: vector for x at each grid point
x = sym('x',[n_grid-1,1]);

%% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms
syms a1 a2 a3 a4 a5 a6 a7 a8 a9

p = [a1,a2,a3,a4,a5,a6,a7,a8,a9];

%% CONSTANTS ( for these no sensitivities will be computed )
% this part is optional and can be ommited

%% SYSTEM EQUATIONS

% create symbolic variable for time
syms t

% time vector
% t_spline = [0, 7.6875, 15.375, 23.0625, 30.75, 38.4375, 46.125, 53.8125, 61.5];

% time-dependent spline for the growth rate
ab = log(am_spline_pos(t,9,0,a1,7.6875,a2,15.375,a3,23.0625,...
    a4,30.75,a5,38.4375,a6,46.125,a7,53.8125,a8,61.5,a9,0,0.0));

xdot = sym(zeros(size(x)));

%% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

%% OBSERVALES
y = ab;

%% SYSTEM STRUCT

model.sym.x = x;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.x0 = x0;
model.sym.y = y;
end