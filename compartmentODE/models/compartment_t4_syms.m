function [model] = compartment_t4_syms( )
% compartment model with time-dependent rates

model.param = 'log';
model.forward = true;
model.adjoint = false;

%% STATES
% create state syms: two compartments
x = sym('x',[2,1]);

%% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms
syms t11 t12 t13 t21 t22 t23 a11 a12 a13 a21 a22 a23

p = [t11,t12,t13,t21,t22,t23,a11,a12,a13,a21,a22,a23];

%% CONSTANTS ( for these no sensitivities will be computed )
% this part is optional and can be ommited

% create parameter syms
k = sym('k',[2,1]);

%% SYSTEM EQUATIONS

% create symbolic variable for time
syms t

% compute splines through parameter dependent nodes
a1t = log(am_spline_pos(t,3,0,a11,30,a12,46.5,a13,0,0.0));
a2t = log(am_spline_pos(t,3,0,a21,30,a22,46.5,a23,0,0.0));
t1t = am_spline_pos(t,3,0,t11,30,t12,46.5,t13,0,0.0);
t2t = am_spline_pos(t,3,0,t21,30,t22,46.5,t23,0,0.0);

xdot = sym(zeros(size(x)));
% A -> B t1t; B -> A t2t; growth A a1t; growth B a2t;

xdot(1) = -t1t*x(1) + t2t*x(2) + a1t*x(1);
xdot(2) = t1t*x(1) - t2t*x(2) + a2t*x(2);

%% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

for i = 1:2
    x0(i) = k(i);
end

%% OBSERVALES
N = sum(x);
y = [x;N];


%% SYSTEM STRUCT

model.sym.x = x;
model.sym.k = k;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.x0 = x0;
model.sym.y = y;
end