%% Globals/Givens
clear
load('KFdata_MODIFIED.mat')
p = 6; n = length(t); I = eye(p); rng(100)
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
% L = 0.5;

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

Q = diag([0.05 0.05 0.1 50 50 0.03]);
% Q = diag([0.005 0.005 0.01 0.05 0.05 0.01]);
% Q = Q*diag([0.0001 0.0001 0.75/2 0.001 0.001 0.75/2]);
% Q(1,3) = 0.0017;  Q(3,1) = Q(1,3);
% Q(2,3) = 0.001;  Q(3,2) = Q(2,3);
% Q(2,1) = -0.013;  Q(1,2) = Q(2,1);
% Q(4,6) = 0.0005;  Q(6,4) = Q(4,6);
Q(5,6) = -0.0005;  Q(6,5) = Q(5,6);
% Q(5,4) = -0.00005;  Q(4,5) = Q(5,4);
Q(3,6) = -0.03;  Q(6,3) = Q(3,6);
% Q(2,5) = 0.000021;  Q(5,2) = Q(3,6);
% Q = zeros(6);
Qkm1 = Q;

% P0 = diag([2.5 2.5 .2 10 10 .1]);
P0 = diag([0.001 0.001 0.01 0.001 0.001 0.01]);


epyk = zeros(1,n-1);
% exk  = zeros(7,n-1);
exk  = zeros(1,n-1);
nis  = zeros(1,n-1);
% nees = zeros(7,n-1);
nees = zeros(1,n-1);

NLdynrecord = zeros(p,n);
NN = 1;

NEESsamps = zeros(NN,n-1);
NISsamps = zeros(NN,n-1);

xkrecord = zeros(p,n);
nldyn = zeros(p,n);
pk = zeros(p,n);
pkrecord = zeros(p,n);
EXKrec = zeros(p,n-1);
exkrec = zeros(p,n-1);
EYKrec = zeros(p-1,n-1);
eykrec = zeros(p-1,n-1);
sk = zeros(p-1,n-1);
skrecord = zeros(p-1,n-1);
for test=1:NN
    
    xtrue = zeros(p,n);
%     xtrue(:,1) = x0;
    xtrue(:,1) = mvnrnd(x0,P0);
    
    for k=2:n
        
        W = mvnrnd(zeros(6,1),Q)';
%         W = zeros(6,1);
        xtrue(:,k) = correct(NLdyn(xtrue(:,k-1),u,W));
%         wrec(:,k) = W;
    end
    
    xk = zeros(p,n);
%     xk(:,1) = mvnrnd(x0,P0);
    xk(:,1) = x0;
    Ppkm1 = diag([100 100 2 1000 1000 2]);
    
    for k=2:n
        % Lookup full state to linearize about
        xig = xk(1,k-1); etag = xk(2,k-1); thetag = xk(3,k-1);
        xia = xk(4,k-1); etaa = xk(5,k-1); thetaa = xk(6,k-1);
        
        % Jacobian estimate
        Fkm1 = Fpred(vg,thetag,va,thetaa);
        
        % Estimation-Error Covariance
        Pmk = Fkm1*Ppkm1*Fkm1' + Qkm1;
        
        
        
        % State Estimate
%         W = mvnrnd(zeros(6,1),Q)';
        xmk = correct(NLdyn(xk(:,k-1),u,zeros(6,1)));
        % record NL state estimate
        nldyn(:,k) = xmk;
        
        % Jacobian estimate
        Hk = Hpred(xmk(1),xmk(2),xmk(4),xmk(5));
        
        % Kalman Gain
        Kk  = Pmk*Hk' / (Hk*Pmk*Hk' + R);
        
        % Nonlinear Measurement Innovation
        ymkp1 = NLmeas(xmk);
        ymkp1(1) = wrapToPi(ymkp1(1)); ymkp1(3) = wrapToPi(ymkp1(3));
%         PRED(:,k) = ymkp1;
%         eyk = [wrapToPi(ydata(1,k)-ymkp1(1)); ...
%             ydata(2,k)-ymkp1(2); ...
%             wrapToPi(ydata(3,k)-ymkp1(3)); ...
%             ydata(4:5,k)-ymkp1(4:5)];
        eyk = [-angdiff(ydata(1,k),ymkp1(1)); ...
            ydata(2,k)-ymkp1(2); ...
            -angdiff(ydata(3,k),ymkp1(3)); ...
            ydata(4:5,k)-ymkp1(4:5)];
        eykrec(:,k-1) = eyk;
        
        % NIS
        epyk(k-1) = NIS(eyk,Hk,Pmk,R);
        
        % Update State Estimate
        xk(:,k) = correct(xmk + Kk*eyk);
        
        % Update Estimation-Error Covariance
        Ppk = (I - Kk*Hk)*Pmk;
        pk(:,k) = [2*sqrt(Ppk(1,1)) 2*sqrt(Ppk(2,2)) 2*sqrt(Ppk(3,3))...
            2*sqrt(Ppk(4,4)) 2*sqrt(Ppk(5,5)) 2*sqrt(Ppk(6,6))]';
        
        sk(:,k-1) = 2*sqrt(diag(Hk*Pmk*Hk' + R));
        
        
        
        % NEES
        diff = [xtrue(1:2,k)-xk(1:2,k); -angdiff(xtrue(3,k),xk(3,k)); ...
            xtrue(4:5,k)-xk(4:5,k); -angdiff(xtrue(6,k),xk(6,k))];
        exk(k-1) = NEES(diff,Ppk);
        exkrec(:,k-1) = diff;
%         exk(:,k-1) = [NEES(diff(1),Ppk(1,1));...
%             NEES(diff(2),Ppk(2,2));...
%             NEES(diff(3),Ppk(3,3));...
%             NEES(diff(4),Ppk(4,4));...
%             NEES(diff(5),Ppk(5,5));...
%             NEES(diff(6),Ppk(6,6));...
%             NEES(diff,Ppk)];
        Ppkm1 = Ppk;
    end
    EYKrec = EYKrec + eykrec;
    EXKrec = EXKrec + exkrec;
    xkrecord = xkrecord + xk;
    pkrecord = pkrecord + pk;
    skrecord = skrecord + sk;
    NLdynrecord = NLdynrecord + nldyn;
    
    NEESsamps(test,:) = exk;
    NISsamps(test,:) = epyk;
    
    nis = nis+epyk;
    nees = nees+exk;
    if test<NN
        clear xk
    end
end
EYKrec = EYKrec/NN;
EXKrec = EXKrec/NN;
xkrecord = xkrecord/NN;
pkrecord = pkrecord/NN;
skrecord = skrecord/NN;
NLdynrecord = NLdynrecord/NN;
nis = nis/NN;
nees = nees/NN;
NEESS = mean(NEESsamps,1);
NISS = mean(NISsamps,1);
% xk = xkrecord;
% xtrue = NLdynrecord;
% pk = pkrecord;

%% Plot pred angel

% figure
% sgtitle('NL Meas.','fontsize',24,'interpreter','latex')
% 
% subplot(5,1,1); hold on; grid on; grid minor
% plot(t,ydata(1,:),'b','Linewidth',1.5)
% plot(t,PRED(1,:),'Linewidth',2)
% legend('ydata','computed')
% xlabel('Time [s]','fontsize',16,'interpreter','latex')
% ylabel('$\gamma_{ag}$ [rads]','fontsize',20,'interpreter','latex')
% hold off
% 
% subplot(5,1,2); hold on; grid on; grid minor
% plot(t,ydata(2,:),'b','Linewidth',1.5)
% plot(t,PRED(2,:),'Linewidth',2)
% legend('ydata','computed')
% xlabel('Time [s]','fontsize',16,'interpreter','latex')
% ylabel('$\rho_{ga}$ [m]','fontsize',20,'interpreter','latex')
% hold off
% 
% subplot(5,1,3); hold on; grid on; grid minor
% plot(t,ydata(3,:),'b','Linewidth',1.5)
% plot(t,PRED(3,:),'Linewidth',2)
% legend('ydata','computed')
% xlabel('Time [s]','fontsize',16,'interpreter','latex')
% ylabel('$\gamma_{ga}$ [rads]','fontsize',20,'interpreter','latex')
% hold off
% 
% subplot(5,1,4); hold on; grid on; grid minor
% plot(t,ydata(4,:),'b','Linewidth',1.5)
% plot(t,PRED(4,:),'Linewidth',2)
% legend('ydata','computed')
% xlabel('Time [s]','fontsize',16,'interpreter','latex')
% ylabel('$\xi_a$ [m]','fontsize',20,'interpreter','latex')
% hold off
% 
% subplot(5,1,5); hold on; grid on; grid minor
% plot(t,ydata(5,:),'b','Linewidth',1.5)
% plot(t,PRED(5,:),'Linewidth',2)
% legend('ydata','computed')
% xlabel('Time [s]','fontsize',16,'interpreter','latex')
% ylabel('$\eta_a$ [m]','fontsize',20,'interpreter','latex')
% hold off


%% Overhead
% figure; hold on
% xlim([min([xk(1,:),xk(4,:)],[],'all') max([xk(1,:),xk(4,:)],[],'all')])
% ylim([min([xk(2,:),xk(5,:)],[],'all') max([xk(2,:),xk(5,:)],[],'all')])
% for i=3:n
%     plot(xk(1,i-1:i),xk(2,i-1:i),'k','Linewidth',1)
% %     plot(NLdynrecord(1,i-1:i),NLdynrecord(2,i-1:i),'r','Linewidth',0.5)
%     
%     plot(xk(4,i-1:i),xk(5,i-1:i),'b','Linewidth',2)
% %     plot(NLdynrecord(4,i-1:i),NLdynrecord(5,i-1:i),'r','Linewidth',0.5)
%     pause(.001)
% end

%% EXK Errors

figure
sgtitle('$e_{x,k}$ with $2\sigma$ bounds from $P_k^+$','fontsize',24,'interpreter','latex')

subplot(6,1,1); hold on; grid on; grid minor
plot(t(2:end),EXKrec(1,:),'Linewidth',2)
plot(t(35:end),-pkrecord(1,35:end),'k--','Linewidth',2)
plot(t(35:end),pkrecord(1,35:end),'k--','Linewidth',2)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\xi_g$ error [m]','fontsize',20,'interpreter','latex')
hold off

subplot(6,1,2); hold on; grid on; grid minor
plot(t(2:end),EXKrec(2,:),'Linewidth',2)
plot(t(35:end),-pkrecord(2,35:end),'k--','Linewidth',2)
plot(t(35:end),pkrecord(2,35:end),'k--','Linewidth',2)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\eta_g$ error [m]','fontsize',20,'interpreter','latex')
hold off

subplot(6,1,3); hold on; grid on; grid minor
plot(t(2:end),EXKrec(3,:),'Linewidth',2)
plot(t(2:end),-pkrecord(3,2:end),'k--','Linewidth',2)
plot(t(2:end),pkrecord(3,2:end),'k--','Linewidth',2)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\theta_g$ error [rad]','fontsize',20,'interpreter','latex')
hold off

subplot(6,1,4); hold on; grid on; grid minor
plot(t(2:end),EXKrec(4,:),'Linewidth',2)
plot(t(2:end),-pkrecord(4,2:end),'k--','Linewidth',2)
plot(t(2:end),pkrecord(4,2:end),'k--','Linewidth',2)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\xi_a$ error [m]','fontsize',20,'interpreter','latex')
hold off

subplot(6,1,5); hold on; grid on; grid minor
plot(t(2:end),EXKrec(5,:),'Linewidth',2)
plot(t(2:end),-pkrecord(5,2:end),'k--','Linewidth',2)
plot(t(2:end),pkrecord(5,2:end),'k--','Linewidth',2)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\eta_a$ error [m]','fontsize',20,'interpreter','latex')
hold off

subplot(6,1,6); hold on; grid on; grid minor
plot(t(2:end),EXKrec(6,:),'Linewidth',2)
plot(t(2:end),-pkrecord(6,2:end),'k--','Linewidth',2)
plot(t(2:end),pkrecord(6,2:end),'k--','Linewidth',2)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\theta_a$ error [rad]','fontsize',20,'interpreter','latex')
hold off


%% EYK Errors

figure
sgtitle('$e_{y,k}$ with $2\sigma$ bounds from $S_k^+$','fontsize',24,'interpreter','latex')

subplot(5,1,1); hold on; grid on; grid minor
plot(t(2:end),EYKrec(1,:),'Linewidth',2)
plot(t(2:end),-skrecord(1,:),'k--','Linewidth',2)
plot(t(2:end),skrecord(1,:),'k--','Linewidth',2)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\gamma_{ag}$ [m]','fontsize',20,'interpreter','latex')
hold off

subplot(5,1,2); hold on; grid on; grid minor
plot(t(2:end),EYKrec(2,:),'Linewidth',2)
plot(t(2:end),-skrecord(2,:),'k--','Linewidth',2)
plot(t(2:end),skrecord(2,:),'k--','Linewidth',2)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\rho_{ga}$ [m]','fontsize',20,'interpreter','latex')
hold off

subplot(5,1,3); hold on; grid on; grid minor
plot(t(2:end),EYKrec(3,:),'Linewidth',2)
plot(t(2:end),-skrecord(3,:),'k--','Linewidth',2)
plot(t(2:end),skrecord(3,:),'k--','Linewidth',2)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\gamma_{ga}$ [m]','fontsize',20,'interpreter','latex')
hold off

subplot(5,1,4); hold on; grid on; grid minor
plot(t(2:end),EYKrec(4,:),'Linewidth',2)
plot(t(2:end),-skrecord(4,:),'k--','Linewidth',2)
plot(t(2:end),skrecord(4,:),'k--','Linewidth',2)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\xi_a$ [m]','fontsize',20,'interpreter','latex')
hold off

subplot(5,1,5); hold on; grid on; grid minor
plot(t(2:end),EYKrec(5,:),'Linewidth',2)
plot(t(2:end),-skrecord(5,:),'k--','Linewidth',2)
plot(t(2:end),skrecord(5,:),'k--','Linewidth',2)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\eta_a$ [m]','fontsize',20,'interpreter','latex')
hold off


%% EKF States
% States vs Time, Full NL Dynamics Simulation
% t = t(40:end);
% xk = xk(:,40:end);
% ydata = ydata(:,40:end);
% xk = xkrecord;
% xtrue = NLdynrecord;

figure
sgtitle('Simulated Ground Truth','fontsize',24,'interpreter','latex')

subplot(6,1,1); hold on; grid on; grid minor
plot(t,xtrue(1,:),'b','Linewidth',1)
plot(t,xk(1,:),'Linewidth',2)
plot(t(10:end),xk(1,10:end)-pk(1,10:end),'b--','Linewidth',1)
plot(t(10:end),xk(1,10:end)+pk(1,10:end),'b--','Linewidth',1)
% legend('NLdyn','filter')
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\xi_g$ [m]','fontsize',20,'interpreter','latex')
hold off

subplot(6,1,2); hold on; grid on; grid minor
plot(t,xtrue(2,:),'b','Linewidth',1)
plot(t,xk(2,:),'Linewidth',2)
plot(t(10:end),xk(2,10:end)-pk(2,10:end),'b--','Linewidth',1)
plot(t(10:end),xk(2,10:end)+pk(2,10:end),'b--','Linewidth',1)
% legend('NLdyn','filter')
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\eta_g$ [m]','fontsize',20,'interpreter','latex')
hold off

subplot(6,1,3); hold on; grid on; grid minor
plot(t,xtrue(3,:),'b.','Linewidth',1)
plot(t,xk(3,:),'Linewidth',2)
plot(t(10:end),xk(3,10:end)-pk(3,10:end),'b--','Linewidth',1)
plot(t(10:end),xk(3,10:end)+pk(3,10:end),'b--','Linewidth',1)
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\theta_g$ [rad]','fontsize',20,'interpreter','latex')
hold off

subplot(6,1,4); hold on; grid on; grid minor
plot(t,xtrue(4,:),'b','Linewidth',1)
plot(t,xk(4,:),'Linewidth',2)
plot(t(10:end),xk(4,10:end)-pk(4,10:end),'b--','Linewidth',1)
plot(t(10:end),xk(4,10:end)+pk(4,10:end),'b--','Linewidth',1)
% plot(t,ydata(4,:))
% legend('NLdyn','filter','meas.')
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\xi_a$ [m]','fontsize',20,'interpreter','latex')
hold off

subplot(6,1,5); hold on; grid on; grid minor
plot(t,xtrue(5,:),'b','Linewidth',1)
plot(t,xk(5,:),'Linewidth',2)
plot(t(10:end),xk(5,10:end)-pk(5,10:end),'b--','Linewidth',1)
plot(t(10:end),xk(5,10:end)+pk(5,10:end),'b--','Linewidth',1)
% plot(t,ydata(5,:))
% legend('NLdyn','filter','meas.')
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\eta_a$ [m]','fontsize',20,'interpreter','latex')
hold off

subplot(6,1,6); hold on; grid on; grid minor
plot(t,xtrue(6,:),'b.','Linewidth',1)
plot(t,xk(6,:),'Linewidth',2)
plot(t(10:end),xk(6,10:end)-pk(6,10:end),'b--','Linewidth',1)
plot(t(10:end),xk(6,10:end)+pk(6,10:end),'b--','Linewidth',1)
% legend('NLdyn','filter')
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\theta_a$ [rad]','fontsize',20,'interpreter','latex')
hold off

%% Noisy Measurements
nlmeas = zeros(5,n);

for k=1:length(xk)
    nlmeas(:,k) = NLmeas(xk(:,k));
    nlmeas(1,k) = wrapToPi(nlmeas(1,k));
    nlmeas(3,k) = wrapToPi(nlmeas(3,k));
end


figure
sgtitle('Simulated Data','fontsize',24,'interpreter','latex')

subplot(5,1,1); hold on; grid on; grid minor
plot(t,ydata(1,:),'b','Linewidth',1.5)
plot(t,nlmeas(1,:),'Linewidth',2)
% legend('ydata','computed')
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\gamma_{ag}$ [rads]','fontsize',20,'interpreter','latex')
hold off

subplot(5,1,2); hold on; grid on; grid minor
plot(t,ydata(2,:),'b','Linewidth',1.5)
plot(t,nlmeas(2,:),'Linewidth',2)
% legend('ydata','computed')
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\rho_{ga}$ [m]','fontsize',20,'interpreter','latex')
hold off

subplot(5,1,3); hold on; grid on; grid minor
plot(t,ydata(3,:),'b','Linewidth',1.5)
plot(t,nlmeas(3,:),'Linewidth',2)
% legend('ydata','computed')
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\gamma_{ga}$ [rads]','fontsize',20,'interpreter','latex')
hold off

subplot(5,1,4); hold on; grid on; grid minor
plot(t,ydata(4,:),'b','Linewidth',1.5)
plot(t,nlmeas(4,:),'Linewidth',2)
% legend('ydata','computed')
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\xi_a$ [m]','fontsize',20,'interpreter','latex')
hold off

subplot(5,1,5); hold on; grid on; grid minor
plot(t,ydata(5,:),'b','Linewidth',1.5)
plot(t,nlmeas(5,:),'Linewidth',2)
% legend('ydata','computed')
xlabel('Time [s]','fontsize',16,'interpreter','latex')
ylabel('$\eta_a$ [m]','fontsize',20,'interpreter','latex')
hold off

%% NEES

alpha = 0.05;
r1 = chi2inv(alpha/2,NN*6)/NN;
r2 = chi2inv(1- alpha/2,NN*6)/NN;

figure;
subplot(2,1,1)
hold on
plot(1:length(nees),nees,'.')
plot(1:length(nees),r1*ones(1,length(nees)),'r--','Linewidth',2)
plot(1:length(nees),r2*ones(1,length(nees)),'r--','Linewidth',2)
xlabel('Time step, k','fontsize',16,'interpreter','latex')
ylabel('$\bar{\epsilon}_x$ [m]','fontsize',20,'interpreter','latex')
title('NEES Estimation Results','fontsize',20,'interpreter','latex')
% title(diag(Q)')
hold off


% r1 = chi2inv(alpha/2,NN)/NN;
% r2 = chi2inv(1- alpha/2,NN)/NN;
% 
% 
% figure;
% subplot(6,1,1); hold on
% plot(1:n-1,nees(1,:),'.')
% plot(1:n-1,r1*ones(1,n-1),'r--')
% plot(1:n-1,r2*ones(1,n-1),'r--')
% xlabel('Time step, k','fontsize',14,'interpreter','latex')
% ylabel('$\bar{\epsilon}_x$ [m]','fontsize',28,'interpreter','latex')
% % title(diag(Q)')
% hold off
% 
% subplot(6,1,2); hold on
% plot(1:n-1,nees(2,:),'.')
% plot(1:n-1,r1*ones(1,n-1),'r--')
% plot(1:n-1,r2*ones(1,n-1),'r--')
% xlabel('Time step, k','fontsize',14,'interpreter','latex')
% ylabel('$\bar{\epsilon}_x$ [m]','fontsize',28,'interpreter','latex')
% % title(diag(Q)')
% hold off
% 
% subplot(6,1,3); hold on
% plot(1:n-1,nees(3,:),'.')
% plot(1:n-1,r1*ones(1,n-1),'r--')
% plot(1:n-1,r2*ones(1,n-1),'r--')
% xlabel('Time step, k','fontsize',14,'interpreter','latex')
% ylabel('$\bar{\epsilon}_x$ [m]','fontsize',28,'interpreter','latex')
% % title(diag(Q)')
% hold off
% 
% subplot(6,1,4); hold on
% plot(1:n-1,nees(4,:),'.')
% plot(1:n-1,r1*ones(1,n-1),'r--')
% plot(1:n-1,r2*ones(1,n-1),'r--')
% xlabel('Time step, k','fontsize',14,'interpreter','latex')
% ylabel('$\bar{\epsilon}_x$ [m]','fontsize',28,'interpreter','latex')
% % title(diag(Q)')
% hold off
% 
% subplot(6,1,5); hold on
% plot(1:n-1,nees(5,:),'.')
% plot(1:n-1,r1*ones(1,n-1),'r--')
% plot(1:n-1,r2*ones(1,n-1),'r--')
% xlabel('Time step, k','fontsize',14,'interpreter','latex')
% ylabel('$\bar{\epsilon}_x$ [m]','fontsize',28,'interpreter','latex')
% % title(diag(Q)')
% hold off
% 
% subplot(6,1,6); hold on
% plot(1:n-1,nees(6,:),'.')
% plot(1:n-1,r1*ones(1,n-1),'r--')
% plot(1:n-1,r2*ones(1,n-1),'r--')
% xlabel('Time step, k','fontsize',14,'interpreter','latex')
% ylabel('$\bar{\epsilon}_x$ [m]','fontsize',28,'interpreter','latex')
% % title(diag(Q)')
% hold off

%% NIS

alpha = 0.05;
r1 = chi2inv(alpha/2,NN*5)/NN;
r2 = chi2inv(1- alpha/2,NN*5)/NN;

subplot(2,1,2)
hold on
plot(1:length(nis),nis,'.')
plot(1:length(nis),r1*ones(1,length(nis)),'r--','Linewidth',2)
plot(1:length(nis),r2*ones(1,length(nis)),'r--','Linewidth',2)
xlabel('Time step, k','fontsize',16,'interpreter','latex')
ylabel('$\bar{\epsilon}_y$ [m]','fontsize',20,'interpreter','latex')
title('NIS Estimation Results','fontsize',20,'interpreter','latex')
hold off
