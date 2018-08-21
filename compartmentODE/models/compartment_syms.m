function [model] = compartment_syms( )
% compartment model with constant rates

model.param = 'log';
model.forward = true;
model.adjoint = false;

%% COMPARTMENT STATES
% create state syms: two compartments
x = sym('x',[2,1]);

%% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms
syms t1 t2 a1 a2

p = [t1,t2,a1,a2];

%% CONSTANTS ( for these no sensitivities will be computed )
% this part is optional and can be ommited

% create parameter syms
k = sym('k',[2,1]);

%% SYSTEM EQUATIONS

% create symbolic variable for time
syms t

xdot = sym(zeros(size(x)));
% ODE for the reactions
% A -> B t1; B -> A t2; growth A a1; growth B a2;

xdot(1) = -t1*x(1) + t2*x(2) + log(a1)*x(1);
xdot(2) = t1*x(1) - t2*x(2) + log(a2)*x(2);

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