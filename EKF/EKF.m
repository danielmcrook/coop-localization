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

Ppkm1 = 1000000000000*ones(p);



Lkm1 = I;
Mk = eye(5);
Qkm1 = 1000*Q; Rk = R;

% xk = zeros(p,n);
% pk = zeros(p,n);
xtrue = zeros(p,n);
xtrue(:,1) = x0;
% xk(:,1) = x0;
epyk = zeros(1,n-1);
exk  = zeros(1,n-1);
nis  = zeros(1,n-1);
nees = zeros(1,n-1);
NLdynrecord = zeros(p,n);
NN = 100;
xkrecord = zeros(p,n);
for test=1:NN
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
        if test==1
            NLdynrecord(:,k) = correct(NLdyn(xk(:,k-1),u));
        end
        xmk = NLdynrecord(:,k);
    %     xmk = correct(xmk);

        % Jacobian estimate
        Hk = Hpred(xmk(1),xmk(2),xmk(4),xmk(5));

        % Kalman Gain
        Kk  = Pmk*Hk' / (Hk*Pmk*Hk' + Mk*Rk*Mk');

        % Nonlinear Measurement Innovation
        eykp1  = ydata(:,k) - NLmeas(xmk);




        % NIS
        epyk(k-1) = NIS(eykp1,Hk,Pmk,R);

        % Update State Estimate
        xk(:,k) = correct(xmk + correct(Kk*eykp1));
        % Update Estimation-Error Covariance
        Ppk = (I - Kk*Hk)*Pmk;
%         pk(:,k) = [2*sqrt(Ppk(1,1)) 2*sqrt(Ppk(2,2)) 2*sqrt(Ppk(3,3))...
%             2*sqrt(Ppk(4,4)) 2*sqrt(Ppk(5,5)) 2*sqrt(Ppk(6,6))]';

%         xk(:,k) = correct(xk(:,k));
        if test==1
            xtrue(:,k) = correct(NLdyn(xtrue(:,k-1),u));
        end

        % NEES
        exk(k-1) = NEES(xtrue(:,k)-xk(:,k),Ppk);
        
        Ppkm1 = Ppk;
    end
    
    xkrecord = xkrecord + xk;
    nis = nis+epyk;
    nees = nees+exk;
    if test<NN
        clear xk
    end
end
xkrecord = xkrecord/NN;
nis = nis/NN;
nees = nees/NN;
xk = xkrecord;
% for i=1:n
%     xk(:,i) = correct(xk(:,i));
% end


% figure; hold on
% xlim([min([xk(1,:),xk(4,:)],[],'all') max([xk(1,:),xk(4,:)],[],'all')])
% ylim([min([xk(2,:),xk(5,:)],[],'all') max([xk(2,:),xk(5,:)],[],'all')])
% for i=3:n
%     plot(xk(1,i-1:i),xk(2,i-1:i),'k','Linewidth',0.5)
%     plot(NLdynrecord(1,i-1:i),NLdynrecord(2,i-1:i),'r','Linewidth',0.5)
%     
%     plot(xk(4,i-1:i),xk(5,i-1:i),'b','Linewidth',0.5)
%     plot(NLdynrecord(4,i-1:i),NLdynrecord(5,i-1:i),'r','Linewidth',0.5)
%     pause(.001)
% end


%% EKF States
% States vs Time, Full NL Dynamics Simulation
% t = t(40:end);
% xk = xk(:,40:end);
% ydata = ydata(:,40:end);
figure
sgtitle('Filter vs Measurments','fontsize',20,'interpreter','latex')

subplot(6,1,1); hold on; grid on; grid minor
plot(t,xtrue(1,:),'b','Linewidth',1)
plot(t,xk(1,:),'k','Linewidth',1.35)
% plot(t(10:end),xk(1,10:end)-pk(1,10:end),'b--','Linewidth',1)
% plot(t(10:end),xk(1,10:end)+pk(1,10:end),'b--','Linewidth',1)
% ylim([10 20]);
legend('NLdyn','filter')
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\xi_g$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,2); hold on; grid on; grid minor
plot(t,xtrue(2,:),'b','Linewidth',1)
plot(t,xk(2,:),'k','Linewidth',1.35)
% plot(t(10:end),xk(2,10:end)-pk(2,10:end),'b--','Linewidth',1)
% plot(t(10:end),xk(2,10:end)+pk(2,10:end),'b--','Linewidth',1)
% ylim([-5 5]);
legend('NLdyn','filter')
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\eta_g$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,3); hold on; grid on; grid minor
plot(t,xtrue(3,:),'b','Linewidth',1)
plot(t,xk(3,:),'k.','Linewidth',1.35)
% plot(t(10:end),xk(3,10:end)-pk(3,10:end),'b--','Linewidth',1)
% plot(t(10:end),xk(3,10:end)+pk(3,10:end),'b--','Linewidth',1)
% ylim([-5 5]);
legend('NLdyn','filter')
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\theta_g$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,4); hold on; grid on; grid minor
plot(t,xtrue(4,:),'b','Linewidth',1)
plot(t,xk(4,:),'k','Linewidth',1.35)
% plot(t(10:end),xk(4,10:end)-pk(4,10:end),'b--','Linewidth',1)
% plot(t(10:end),xk(4,10:end)+pk(4,10:end),'b--','Linewidth',1)
plot(t,ydata(4,:))
legend('NLdyn','filter','meas.')
% ylim([-200 200]);
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\xi_a$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,5); hold on; grid on; grid minor
plot(t,xtrue(5,:),'b','Linewidth',1)
plot(t,xk(5,:),'k','Linewidth',1.35)
% plot(t(10:end),xk(5,10:end)-pk(5,10:end),'b--','Linewidth',1)
% plot(t(10:end),xk(5,10:end)+pk(5,10:end),'b--','Linewidth',1)
plot(t,ydata(5,:))
legend('NLdyn','filter','meas.')
% ylim([-200 200]);
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\eta_a$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(6,1,6); hold on; grid on; grid minor
plot(t,xtrue(6,:),'b','Linewidth',1)
plot(t,xk(6,:),'k.','Linewidth',1.35)
% plot(t(10:end),xk(6,10:end)-pk(6,10:end),'b--','Linewidth',1)
% plot(t(10:end),xk(6,10:end)+pk(6,10:end),'b--','Linewidth',1)
% ylim([-5 5]);
legend('NLdyn','filter')
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\theta_a$ [m]','fontsize',16,'interpreter','latex')
hold off

%% Noisy Measurements
figure
sgtitle('Noisy Measurements','fontsize',20,'interpreter','latex')

subplot(5,1,1); hold on; grid on; grid minor
plot(t,ydata(1,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\gamma_{ag}$ [rads]','fontsize',16,'interpreter','latex')
hold off

subplot(5,1,2); hold on; grid on; grid minor
plot(t,ydata(2,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\rho_{ga}$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(5,1,3); hold on; grid on; grid minor
plot(t,ydata(3,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\gamma_{ga}$ [rads]','fontsize',16,'interpreter','latex')
hold off

subplot(5,1,4); hold on; grid on; grid minor
plot(t,ydata(4,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\xi_a$ [m]','fontsize',16,'interpreter','latex')
hold off

subplot(5,1,5); hold on; grid on; grid minor
plot(t,ydata(5,:),'Linewidth',1.35)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\eta_a$ [m]','fontsize',16,'interpreter','latex')
hold off

%% consistency


%% NEES

alpha = 0.05;
r1 = chi2inv(alpha/2,NN*n)/NN;
r2 = chi2inv(1- alpha/2,NN*n)/NN;

figure; subplot(2,1,1)
hold on
plot(1:length(nees),nees,'.')
plot(1:length(nees),r1*ones(1,length(nees)),'r--')
plot(1:length(nees),r2*ones(1,length(nees)),'r--')
title('NEES')
hold off

%% NIS

% epyk = NIS(minnov,R);

subplot(2,1,2)
hold on
plot(1:length(nis),nis,'.')
title('NIS')
hold off