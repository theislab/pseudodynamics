function [model] = splineHelp_t_syms( )
% simulation file for the evaluation of the spline for the growth rate 
% at different time points 

model.param = 'log';
model.forward = true;
model.adjoint = false;

%% STATES
% create state syms: two compartments
x = sym('x',[2,1]);

%% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms
syms t1 t2 a11 a12 a13 a21 a22 a23

p = [t1,t2,a11,a12,a13,a21,a22,a23];

%% CONSTANTS ( for these no sensitivities will be computed )
% this part is optional and can be ommited

%% SYSTEM EQUATIONS

% create symbolic variable for time
syms t

% time-dependent spline for the growth rates
a1t = log(am_spline_pos(t,3,0,a11,30,a12,46.5,a13,0,0.0));
a2t = log(am_spline_pos(t,3,0,a21,30,a22,46.5,a23,0,0.0));

xdot = sym(zeros(size(x)));

%% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

%% OBSERVALES
y = [a1t;a2t];


%% SYSTEM STRUCT

model.sym.x = x;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.x0 = x0;
model.sym.y = y;
end