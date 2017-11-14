function [cp, dcpdx] = cspline_00(v, vp, x_k, x_kp1, xp)
% calculate cspline to end point 0 and end slope 0 using hermite cubic
% splines, for given initial slope v_p and value v on interval x_k x_kp1

tp = (xp-x_k)./(x_kp1-x_k);
cp = (2*tp.^3-3*tp.^2+1)*v + (tp.^3-2*tp.^2+tp)*vp*(x_kp1-x_k);
dtdx = 1/(x_kp1-x_k);
dcpdx = (6*tp.^2*dtdx-6*tp*dtdx)*v+(3*tp.^2*dtdx-4*tp*dtdx+dtdx)*vp*(x_kp1-x_k);
end
