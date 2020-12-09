function out = NLdyn(x,u,W)
%NLDYN Summary of this function goes here
%   Detailed explanation goes here

tspan = [0,0.1];

[~,y] = ode45(@(t,y) NLcoop(t,y,u),tspan,x);

out = y(end,:)' + W;

end

