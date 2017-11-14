function [varargout] = simpleSplineSym(x,y,xx)
% simpleSplineSym computes the symbolic expression of a natural cubic spline given by the symbolic variables y on the
% equidistant grid x at the grid points xx
% natural cubic splines
% only equidistant
% no check of inputs
% x vector of positions
% y vector of symbolic variables
% xx interpolation point
% yy symbolic formula dependent on y
yy = sym(size(xx));
dyydx = sym(size(xx));
h = x(2)-x(1);
b = 1/h*(y(2:end)-y(1:end-1));
v = 4*h;
u = 6*(b(2:end)-b(1:end-1));
z = [0;inv(diag(v*ones(length(x)-2,1))+diag(h*ones(length(x)-3,1),-1)+diag(h*ones(length(x)-3,1),1))*u;0];
for i=1:length(xx)
    ind = min(find(x<=xx(i),1,'last'),length(x)-1);
    yy(i) = z(ind+1)/(6*h)*(xx(i)-x(ind))^3+z(ind)/(6*h)*(x(ind+1)-xx(i))^3+...
        (y(ind+1)/h-z(ind+1)/6*h)*(xx(i)-x(ind))+(y(ind)/h-h/6*z(ind))*(x(ind+1)-xx(i));
    if nargout == 2
        dyydx(i)=z(ind+1)/(2*h)*(xx(i)-x(ind))^2-z(ind)/(2*h)*(x(ind+1)-xx(i))^2+b(ind)-h/6*(z(ind+1)-z(ind));
    end
end
if nargout == 1
    varargout{1} = yy;
else
    varargout{1} = yy;
    varargout{2} = dyydx;
end
        
end