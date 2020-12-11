function dx = NLcoop(~,y,u,W)
% ode45 NL equations

% States
thetag = y(3);
thetaa = y(6);

% Controls
vg = u(1);
phig = u(2);
va = u(3);
wa = u(4);

% Extra parameters
L = 0.5;

% Switch cases
if wa > pi/6
    wa = pi/6;
elseif wa < -pi/6
    wa = -pi/6;
end

if va > 20
    va = 20;
elseif va < 10
    va = 10;
end

if phig > 5*pi/12
    phig = 5*pi/12;
elseif phig < -5*pi/12
    phig = -5*pi/12;
end

if vg > 3
    vg = 3;
end

% NL EOM
dx = cast([vg*cos(thetag); vg*sin(thetag); vg*tan(phig)/L;...
    va*cos(thetaa); va*sin(thetaa); wa] + W,'double');

end