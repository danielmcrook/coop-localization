%% Globals/Givens
load('KFdata_MODIFIED.mat')
p = 6; n = length(t); I = eye(p);
%% NL Matrices
A = @(vg,thetag,va,thetaa) [0 0 -vg*sin(thetag) 0 0 0;...
    0 0 vg*cos(thetag) 0 0 0;...
    0 0 0 0 0 0;...
    0 0 0 0 0 -va*sin(thetaa);...
    0 0 0 0 0 va*cos(thetaa);...
    0 0 0 0 0 0];

% B = @(thetag,L,phig,vg,thetaa) [cos(thetag) 0 0 0;...
%     sin(thetag) 0 0 0;...
%     tan(phig)/L (vg/L)*sec(phig)^2 0 0;...
%     0 0 cos(thetaa) 0;...
%     0 0 sin(thetaa) 0;...
%     0 0 0 1];

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

% D = zeros(6,4);
%% Nominal Traj. & Linearization
% UGV length, fixed
L = 0.5;

% given initial state
xig = 10;
etag = 0;
thetag = pi/2;
xia = -60;
etaa = 0;
thetaa = -pi/2;

% given controls
vg = 2;
phig = -pi/18;
va = 12;
wa = pi/25;
%% Discretization
DT = 0.1;

Fpred = @(vg,thetag,va,thetaa) I + DT*A(vg,thetag,va,thetaa);
Hpred = @(xig,etag,xia,etaa) C(xig,etag,xia,etaa);

x0 = [xig etag thetag xia etaa thetaa]';
u = [vg phig va wa]';
%% Filter

xpkm1 = x0;
Ppkm1 = 1000000000000*ones(p);



Lkm1 = I;
Mk = eye(5);
Qkm1 = Q; Rk = R;

xk = zeros(p,n);
xk(:,1) = x0;

for k=2:n
    % Lookup full state to linearize about
    xig = xk(1,k-1); etag = xk(2,k-1); thetag = xk(3,k-1);
    xia = xk(4,k-1); etaa = xk(5,k-1); thetaa = xk(6,k-1);
    
    % Jacobian estimate
    Fkm1 = Fpred(vg,thetag,va,thetaa);
    
    
    % Estimation-Error Covariance
    Pmk = Fkm1*Ppkm1*Fkm1' + Lkm1*Qkm1*Lkm1';
    % State Estimate
    xmk = NLdyn(xpkm1,u);
    
    % Jacobian estimate
    Hk = Hpred(xmk(1),xmk(2),xmk(4),xmk(5));
    
    % Kalman Gain
    Kk  = Pmk*Hk' * inv(Hk*Pmk*Hk' + Mk*Rk*Mk');
    
    % Nonlinear Measurement Innovation
    eykp1  = ydata(:,k) - NLmeas(xmk);
    % Update State Estimate
    xk(:,k) = xmk + Kk*eykp1;
    % Update Estimation-Error Covariance
    Ppk = (I - Kk*Hk)*Pmk;
    
    Ppkm1 = Ppk;
    xpkm1 = xk(:,k);
    
    if k==234
        b=1;
    end
    
end

% figure; hold on
% xlim([min([xk(1,:),xk(4,:)],[],'all') max([xk(1,:),xk(4,:)],[],'all')])
% ylim([min([xk(2,:),xk(5,:)],[],'all') max([xk(2,:),xk(5,:)],[],'all')])
% for i=2:n
%     plot(xk(1,i-1:i),xk(2,i-1:i),'k','Linewidth',1.35)
%     plot(xk(4,i-1:i),xk(5,i-1:i),'b','Linewidth',1)
%     pause(.001)
% end


%% plotting
% States vs Time, Full NL Dynamics Simulation
% t = t(40:end);
% xk = xk(:,40:end);
% ydata = ydata(:,40:end);
figure
sgtitle('Filter vs Measurments','fontsize',20,'interpreter','latex')

subplot(6,1,1); hold on; grid on; grid minor
plot(t,xk(1,:),'Linewidth',1.35)
% ylim([10 20]);
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\xi_g$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,2); hold on; grid on; grid minor
plot(t,xk(2,:),'Linewidth',1.35)
% ylim([-5 5]);
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\eta_g$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,3); hold on; grid on; grid minor
plot(t,xk(3,:),'Linewidth',1.35)
% ylim([-5 5]);
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\theta_g$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,4); hold on; grid on; grid minor
plot(t,xk(4,:),'Linewidth',1.35)
plot(t,ydata(4,:))
legend('filter','meas.')
% ylim([-200 200]);
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\xi_a$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,5); hold on; grid on; grid minor
plot(t,xk(5,:),'Linewidth',1.35)
plot(t,ydata(5,:))
legend('filter','meas.')
% ylim([-200 200]);
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\eta_a$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,6); hold on; grid on; grid minor
plot(t,xk(6,:),'Linewidth',1.35)
% ylim([-5 5]);
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\theta_a$ [m]','fontsize',16,'interpreter','latex')
hold off