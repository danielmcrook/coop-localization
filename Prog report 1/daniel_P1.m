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

%% Full Nonlinear Perturbation Dynamics
% tspan = [0,100];
tspan = 0:0.1:100;
x0 = x0 + [0 1 0 0 0 0.1]';
[t,y] = ode45(@(t,y) NLcoop(t,y,u,L),tspan,x0);

% correct rad data
for i=1:length(t)
    thetag = y(i,3);
    thetaa = y(i,6);
    if thetag>pi
        y(i:end,3)=thetag-2*pi;
    elseif thetag<-pi
        y(i:end,3)=y(i:end,3)+2*pi;
    end
    
    if thetaa>pi
        y(i:end,6)=y(i:end,6)-2*pi;
    elseif thetaa<-pi
        y(i:end,6)=thetaa+2*pi;
    end
end

% States vs Time, Full NL Dynamics Simulation
figure
sgtitle('Full Nonlinear Model States','fontsize',20,'interpreter','latex')

subplot(6,1,1); hold on; grid on; grid minor
plot(t,y(:,1),'Linewidth',1.35)
ylim([10 20]);
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\xi_g$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,2); hold on; grid on; grid minor
plot(t,y(:,2),'Linewidth',1.35)
ylim([-5 5]);
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\eta_g$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,3); hold on; grid on; grid minor
plot(t,y(:,3),'Linewidth',1.35)
ylim([-5 5]);
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\theta_g$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,4); hold on; grid on; grid minor
plot(t,y(:,4),'Linewidth',1.35)
ylim([-200 200]);
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\xi_a$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,5); hold on; grid on; grid minor
plot(t,y(:,5),'Linewidth',1.35)
ylim([-200 200]);
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\eta_a$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,6); hold on; grid on; grid minor
plot(t,y(:,6),'Linewidth',1.35)
ylim([-5 5]);
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\theta_a$ [m]','fontsize',16,'interpreter','latex')
hold off

%% Full Nonlinear Model Data Simulation
yk = @(xig,etag,thetag,xia,etaa,thetaa) [atan2((etaa-etag),(xia-xig))-thetag...
    sqrt((xig-xia)^2 + (etag-etaa)^2)...
    atan2((etag-etaa),(xig-xia))-thetaa...
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
    
    % correct rad data
    gammaag = yt(1,i);
    gammaga = yt(3,i);
    if gammaag>pi
        yt(1,i)=gammaag-2*pi;
    elseif gammaag<-pi
        yt(1,i)=gammaag+2*pi;
    end
    
    if gammaga>pi
        yt(3,i)=gammaga-2*pi;
    elseif gammaga<-pi
        yt(3,i)=gammaga+2*pi;
    end
end

yt = yt';

figure
sgtitle('Full Nonlinear Model Measurements','fontsize',20,'interpreter','latex')

subplot(5,1,1); hold on; grid on; grid minor
plot(t,yt(:,1),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\gamma_{ag}$ [rads]','fontsize',16,'interpreter','latex')
hold off

subplot(5,1,2); hold on; grid on; grid minor
plot(t,yt(:,2),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\rho_{ga}$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(5,1,3); hold on; grid on; grid minor
plot(t,yt(:,3),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\gamma_{ga}$ [rads]','fontsize',16,'interpreter','latex')
hold off

subplot(5,1,4); hold on; grid on; grid minor
plot(t,yt(:,4),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\xi_a$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(5,1,5); hold on; grid on; grid minor
plot(t,yt(:,5),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\eta_a$ [m]','fontsize',16,'interpreter','latex')
hold off

%% DT LTI Simulations
x0 = [xig etag thetag xia etaa thetaa]';
[t,y] = ode45(@(t,y) NLcoop(t,y,u,L),tspan,x0);

clear yk
xk = [0 1 0 0 0 0.1]';

vg = 2;
va = 12;
phig = -pi/18;
wa = pi/25;

for k=1:1000
    % Lookup full state to linearize about
    xig = y(k,1);
    etag = y(k,2);
    thetag = y(k,3);
    xia = y(k,4);
    etaa = y(k,5);
    thetaa = y(k,6);
    % Re-Linearize Matrices
    Alin = A(vg,thetag,va,thetaa);
    Clin = C(xig,etag,xia,etaa);
    
    F = eye(6) + DT*Alin;
    H = DT*Clin;
    % DT Dynamical Simulation
    xk(:,k+1) = F*xk(:,k);
    yk(:,k) = H*xk(:,k+1);
end

tk = linspace(0,DT*k,k+1);

% Add perturbations to full state
xL = xk + y';

%% DT LTI
figure
sgtitle('Linearized Approximate States','fontsize',20,'interpreter','latex')

subplot(6,1,1); hold on; grid on; grid minor
plot(tk,xL(1,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\xi_g$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,2); hold on; grid on; grid minor
plot(tk,xL(2,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\eta_g$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,3); hold on; grid on; grid minor
plot(tk,xL(3,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\theta_g$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,4); hold on; grid on; grid minor
plot(tk,xL(4,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\xi_a$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,5); hold on; grid on; grid minor
plot(tk,xL(5,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\eta_a$ [m]','fontsize',16,'interpreter','latex')

subplot(6,1,6); hold on; grid on; grid minor
plot(tk,xL(6,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\theta_a$ [m]','fontsize',16,'interpreter','latex')
hold off

%% DT LTV Perturbations
figure
sgtitle('Linearized Approximate Perturbations','fontsize',20,'interpreter','latex')

subplot(6,1,1); hold on; grid on; grid minor
plot(tk,xk(1,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\delta\xi_g$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,2); hold on; grid on; grid minor
plot(tk,xk(2,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\delta\eta_g$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,3); hold on; grid on; grid minor
plot(tk,xk(3,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\delta\theta_g$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,4); hold on; grid on; grid minor
plot(tk,xk(4,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\delta\xi_a$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,5); hold on; grid on; grid minor
plot(tk,xk(5,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\delta\eta_a$ [m]','fontsize',16,'interpreter','latex')

subplot(6,1,6); hold on; grid on; grid minor
plot(tk,xk(6,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\delta\theta_a$ [m]','fontsize',16,'interpreter','latex')
hold off

%% DT LTI Measurements
tk = tk(2:end);
% Add perturbations to full measurments
yL = yk + yt(2:end,:)';

figure
sgtitle('Linearized Approximate Measurements','fontsize',20,'interpreter','latex')

subplot(5,1,1); hold on; grid on; grid minor
plot(tk,yL(1,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\gamma_{ag}$ [rads]','fontsize',16,'interpreter','latex')
hold off

subplot(5,1,2); hold on; grid on; grid minor
plot(tk,yL(2,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\rho_{ga}$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(5,1,3); hold on; grid on; grid minor
plot(tk,yL(3,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\gamma_{ga}$ [rads]','fontsize',16,'interpreter','latex')
hold off

subplot(5,1,4); hold on; grid on; grid minor
plot(tk,yL(4,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\xi_a$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(5,1,5); hold on; grid on; grid minor
plot(tk,yL(5,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\eta_a$ [m]','fontsize',16,'interpreter','latex')
hold off