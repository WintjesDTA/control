function [xlnotcovered, xrnotcovered] = cover( goalkeeper, attacker, theta)
%COVER Summary of this function goes here
%   Detailed explanation goes here
close all;
goalsize = 7.32;
d = .9;
xk = goalkeeper;
xa = attacker;
xarm1 = xk(1) - cos(theta)*d;
yarm1 = xk(2) - sin(theta)*d;
xarm2 = xk(1) + cos(theta)*d;
yarm2 = xk(2) + sin(theta)*d;

deltax1 = (xa(1)-xarm1)*(xa(2)/(xa(2)-yarm1));
deltax2 = (xa(1)-xarm2)*(xa(2)/(xa(2)-yarm2));
xr = [xa(1) - deltax1, 0];
xl = [xa(1) - deltax2, 0];

if xr(1) < xl(1)
    b = xl;
    xl = xr;
    xr = b;
end
xlnotcovered = abs(min(0,-goalsize/2 - xl(1)));
xrnotcovered = max(0,goalsize/2 - xr(1));
end
