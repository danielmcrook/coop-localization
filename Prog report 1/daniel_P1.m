%% NL Matrix Definitions
A = @(vg,thetag,va,thetaa) [0 0 -vg*sin(thetag) 0 0 0;...
    0 0 vg*cos(thetag) 0 0 0;...
    0 0 0 0 0 0;...
    0 0 0 0 0 -va*sin(thetaa);...
    0 0 0 0 0 va*cos(thetaa);...
    0 0 0 0 0 0];

B = @(thetag,L,phig,vg,thetaa) [cos(thetag) 0 0 0;...
    sin(thetag) 0 0 0;...
    tan(phig)/L (vg/L)*sec(phig)^2 0 0;...
    0 0 cos(thetaa) 0;...
    0 0 sin(thetaa) 0;...
    0 0 0 1];

C = @(xig,etag,xia,etaa) [((etaa-etag)/(xia-xig)^2)/(1+((etaa-etag)/(xia-xig))^2)...
    -(1/(xia-xig))/(1+((etaa-etag)/(xia-xig))^2) -1 ...
    -((etaa-etag)/(xia-xig)^2)/(1+((etaa-etag)/(xia-xig))^2)...
    (1/(xia-xig))/(1+((etaa-etag)/(xia-xig))^2) 0;...
    (xig-xia)/sqrt((xig-xia)^2+(etag-etaa)^2)...
    (etag-etaa)/sqrt((xig-xia)^2+(etag-etaa)^2) 0 ...
    -(xig-xia)/sqrt((xig-xia)^2+(etag-etaa)^2)...
    (etag-etaa)/sqrt((xig-xia)^2+(etag-etaa)^2) 0;...
    -((etag-etaa)/(xig-xia)^2)/(1+((etag-etaa)/(xig-xia))^2)...
    (1/(xig-xia))/(1+((etag-etaa)/(xig-xia))^2) 0 ...
    ((etag-etaa)/(xig-xia)^2)/(1+((etag-etaa)/(xig-xia))^2)...
    -(1/(xig-xia))/(1+((etag-etaa)/(xig-xia))^2) -1;...
    0 0 0 1 0 0; 0 0 0 0 1 0];

D = zeros(6,4);

%% Linearization
L = 0.5;
xig = 10;
etag = 0;
thetag = pi/2;
vg = 2;
phig = -pi/18;
xia = -60;
etaa = 0;
thetaa = -pi/2;
va = 12;
wa = pi/25;

Alin = A(vg,thetag,va,thetaa);
Blin = B(thetag,L,phig,vg,thetaa);
Clin = C(xig,etag,xia,etaa);
Dlin = D;

%% Discretization
DT = 0.1;

z = expm(DT*[Alin Blin; zeros(4,10)]);
F = z(1:6,1:6);
G = z(1:6,7:10);
H = Clin;
M = Dlin;

x0 = [xig etag thetag xia etaa thetaa]';
x = x0;
u = [vg phig va wa]';

%% Observability
O = [H; H*F; H*F^2; H*F^3; H*F^4; H*F^5; H*F^6];
rank(O)

%% Full Nonlinear Perturbation Dynamics
% tspan = [0,100];
tspan = 0:0.01:100;
x0 = x0 + [0 1 0 0 0 0.1]';
[t,y] = ode45(@(t,y) NLcoop(t,y,u,L),tspan,x0);

% correct rad data
for i=1:length(t)
    thetag = y(i,3);
    thetaa = y(i,6);
    if thetag>2*pi
        y(i:end,3)=thetag-2*pi;
    elseif thetag<-2*pi
        y(i:end,3)=y(i:end,3)+2*pi;
    end
    
    if thetaa>2*pi
        y(i:end,6)=y(i:end,6)-2*pi;
    elseif thetaa<-2*pi
        y(i:end,6)=thetaa+2*pi;
    end
end

% States vs Time, Full NL Dynamics Simulation
figure
sgtitle('States vs Time, Full NL Dynamics')

subplot(6,1,1); hold on; grid on; grid minor
plot(t,y(:,1))
xlabel('Time [s]')
ylabel('xig [m]')
hold off

subplot(6,1,2); hold on; grid on; grid minor
plot(t,y(:,2))
xlabel('Time [s]')
ylabel('etag [m]')
hold off

subplot(6,1,3); hold on; grid on; grid minor
plot(t,y(:,3))
xlabel('Time [s]')
ylabel('thetag [m]')
hold off

subplot(6,1,4); hold on; grid on; grid minor
plot(t,y(:,4))
xlabel('Time [s]')
ylabel('xia [m]')
hold off

subplot(6,1,5); hold on; grid on; grid minor
plot(t,y(:,5))
xlabel('Time [s]')
ylabel('etaa [m]')
hold off

subplot(6,1,6); hold on; grid on; grid minor
plot(t,y(:,6))
xlabel('Time [s]')
ylabel('thetaa [m]')
hold off

%% Full Nonlinear Model Data Simulation
yk = @(xig,etag,thetag,xia,etaa,thetaa) [atan((etaa-etag)/(xia-xig)) - thetag...
    sqrt((xig-xia)^2 + (etag-etaa)^2)...
    atan((etag-etaa)/(xig-xia)) - thetaa...
    xia etaa]';
clear yt
for i=1:length(t)
    xig = y(i,1);
    etag = y(i,2);
    thetag = y(i,3);
    xia = y(i,4);
    etaa = y(i,5);
    thetaa = y(i,6);
    yt(:,i) = yk(xig,etag,thetag,xia,etaa,thetaa);
end

yt = yt';

figure
sgtitle('Full Nonlinear Measurements')

subplot(5,1,1); hold on; grid on; grid minor
plot(t,yt(:,1))
xlabel('Time [s]')
ylabel('gamma ag [rads]')
hold off

subplot(5,1,2); hold on; grid on; grid minor
plot(t,yt(:,2))
xlabel('Time [s]')
ylabel('rho ga [m]')
hold off

subplot(5,1,3); hold on; grid on; grid minor
plot(t,yt(:,3))
xlabel('Time [s]')
ylabel('gamma ga [rads]')
hold off

subplot(5,1,4); hold on; grid on; grid minor
plot(t,yt(:,4))
xlabel('Time [s]')
ylabel('xia [m]')
hold off

subplot(5,1,5); hold on; grid on; grid minor
plot(t,yt(:,5))
xlabel('Time [s]')
ylabel('etaa [m]')
hold off

%% DT LTI Simulations
clear yk
xk = x0;

for k=1:1000
    xk(:,k+1) = F*xk(:,k);
    yk(:,k) = H*xk(:,k+1);
end

tk = linspace(0,k/10,k+1);

%% DT LTI
figure
sgtitle('States v. Time, Linearized DT')

subplot(6,1,1); hold on; grid on; grid minor
plot(tk,xk(1,:))
xlabel('Time [s]')
ylabel('xig [m]')
hold off

subplot(6,1,2); hold on; grid on; grid minor
plot(tk,xk(2,:))
xlabel('Time [s]')
ylabel('etag [m]')
hold off

subplot(6,1,3); hold on; grid on; grid minor
plot(tk,xk(3,:))
xlabel('Time [s]')
ylabel('thetag [m]')
hold off

subplot(6,1,4); hold on; grid on; grid minor
plot(tk,xk(4,:))
xlabel('Time [s]')
ylabel('xia [m]')
hold off

subplot(6,1,5); hold on; grid on; grid minor
plot(tk,xk(5,:))
xlabel('Time [s]')
ylabel('etaa [m]')

subplot(6,1,6); hold on; grid on; grid minor
plot(tk,xk(6,:))
xlabel('Time [s]')
ylabel('thetaa [m]')
hold off

%% DT LTI Measurements
tk = tk(2:end);
figure
sgtitle('Linearized DT Measurements')

subplot(5,1,1); hold on; grid on; grid minor
plot(tk,yk(1,:))
xlabel('Time [s]')
ylabel('gamma ag [m]')
hold off

subplot(5,1,2); hold on; grid on; grid minor
plot(tk,yk(2,:))
% plot(t,y(:,2))
xlabel('Time [s]')
ylabel('rho ga [m]')
hold off

subplot(5,1,3); hold on; grid on; grid minor
plot(tk,yk(3,:))
% plot(t,y(:,3))
xlabel('Time [s]')
ylabel('gamma ga [m]')
hold off

subplot(5,1,4); hold on; grid on; grid minor
plot(tk,yk(4,:))
% plot(t,y(:,4))
xlabel('Time [s]')
ylabel('xia [m]')
hold off

subplot(5,1,5); hold on; grid on; grid minor
plot(tk,yk(5,:))
% plot(t,y(:,5))
xlabel('Time [s]')
ylabel('etaa [m]')
hold off