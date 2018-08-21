function [model] = splineHelp_t4_syms( )
% simulation file for the evaluation of the spline for the growth and 
% transition rates at different time points 

model.param = 'log';
model.forward = true;
model.adjoint = false;

%% STATES
% create state syms: vector for x at each grid point
x = sym('x',[2,1]);

%% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms
syms t11 t12 t13 t21 t22 t23 a11 a12 a13 a21 a22 a23

p = [t11,t12,t13,t21,t22,t23,a11,a12,a13,a21,a22,a23];

%% CONSTANTS ( for these no sensitivities will be computed )
% this part is optional and can be ommited

%% SYSTEM EQUATIONS

% create symbolic variable for time
syms t

% time-dependent splines for growth and transition rates
a1t = log(am_spline_pos(t,3,0,a11,30,a12,46.5,a13,0,0.0));
a2t = log(am_spline_pos(t,3,0,a21,30,a22,46.5,a23,0,0.0));
t1t = am_spline_pos(t,3,0,t11,30,t12,46.5,t13,0,0.0);
t2t = am_spline_pos(t,3,0,t21,30,t22,46.5,t23,0,0.0);

xdot = sym(zeros(size(x)));

%% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));


%% OBSERVALES
y = [a1t; a2t; t1t; t2t];


%% SYSTEM STRUCT

model.sym.x = x;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.x0 = x0;
model.sym.y = y;
end