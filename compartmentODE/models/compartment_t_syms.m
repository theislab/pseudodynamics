function [model] = compartment_t_syms( )
% compartment model with constant transition and time-dependent growth rates

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

% create parameter syms
k = sym('k',[2,1]);

%% SYSTEM EQUATIONS

% create symbolic variable for time
syms t

% compute splines through parameter dependent nodes
a1t = log(am_spline_pos(t,3,0,a11,30,a12,46.5,a13,0,0.0));
a2t = log(am_spline_pos(t,3,0,a21,30,a22,46.5,a23,0,0.0));

xdot = sym(zeros(size(x)));
% A -> B t1; B -> A t2; growth A a1t; growth B a2t;

xdot(1) = -t1*x(1) + t2*x(2) + a1t*x(1);
xdot(2) = t1*x(1) - t2*x(2) + a2t*x(2);

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