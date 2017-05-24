function [xlnotcovered, xrnotcovered] = coverDraw( goalkeeper, attacker, theta)
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
% Plot soccer field
hold on;
plot([xa(1),xl(1)],[xa(2),xl(2)],'-');
plot([xa(1),xr(1)],[xa(2),xr(2)],'-');
plot([xarm1,xarm2],[yarm1,yarm2]);
leftpole = -7.32/2 - 16.5;
rightpole = -leftpole;
xgoal = [leftpole, leftpole, rightpole, rightpole];
ygoal = [0,16.5,16.5,0];
plot(xgoal,ygoal);
axis([-25 25 0 35]);
plot(goalkeeper(1),goalkeeper(2),'b+');
plot(attacker(1),attacker(2),'b+');
% Goal posts
plot([-goalsize/2,-goalsize/2],[0,3],'black')
plot([goalsize/2,goalsize/2],[0,3],'black')
title('Plot of the lightsoccer situation');
end

