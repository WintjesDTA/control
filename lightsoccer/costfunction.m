function [ J ] = costfunction(keeperx,keepery,theta)
%J Summary of this function goes here
%   Detailed explanation goes here
attacker = [-5,17.5];
[xl, xr] = cover([keeperx,keepery],attacker,theta);
J = xl^2 + xr^2;
end

