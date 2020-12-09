function [Observation,Reward,IsDone,LoggedSignals] = EKFtrain(Action,LoggedSignals)
%% Globals/Givens
load('RLdata.mat','t','R','ydata')
% t = 0:0.1:100;
% R = diag([0.0225 64 0.04 36 36]);
p = 6; n = length(t); I = eye(p); rng(100)
%% NL Matrices
A = @(vg,thetag,va,thetaa) [0 0 -vg*sin(thetag) 0 0 0;...
    0 0 vg*cos(thetag) 0 0 0;...
    0 0 0 0 0 0;...
    0 0 0 0 0 -va*sin(thetaa);...
    0 0 0 0 0 va*cos(thetaa);...
    0 0 0 0 0 0];

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
%% Nominal Traj. & Linearization

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
% Action = [ones(1,6), zeros(1,15), ones(1,6)];
Q = diag(Action(1:6));

Q(1,2) = Action(7);  Q(2,1) = Q(1,2);
Q(1,3) = Action(8);  Q(3,1) = Q(1,3);
Q(1,4) = Action(9);  Q(4,1) = Q(1,4);
Q(1,5) = Action(10);  Q(5,1) = Q(1,5);
Q(1,6) = Action(11);  Q(6,1) = Q(1,6);

Q(2,3) = Action(12);  Q(3,2) = Q(2,3);
Q(2,4) = Action(13);  Q(4,2) = Q(2,4);
Q(2,5) = Action(14);  Q(5,2) = Q(2,5);
Q(2,6) = Action(15);  Q(6,2) = Q(2,6);

Q(3,4) = Action(16);  Q(4,3) = Q(3,4);
Q(3,5) = Action(17);  Q(5,3) = Q(3,5);
Q(3,6) = Action(18);  Q(6,3) = Q(3,6);

Q(4,5) = Action(19);  Q(5,4) = Q(4,5);
Q(4,6) = Action(20);  Q(6,4) = Q(4,6);

Q(5,6) = Action(21);  Q(6,5) = Q(5,6);

[~,P] = cholcov(Q);
% disp(P)
if P > 0
    NN = 0;
    IsDone = 1;
    Reward = -Inf;
    Observation = -100000000000*ones(7,1000);
else
    NN = 10;
end


P0 = diag(Action(22:27));


% epyk = zeros(1,n-1);
% exk  = zeros(1,n-1);
% nis  = zeros(1,n-1);
nees = zeros(1,n-1);

% NLdynrecord = zeros(p,n);


% xkrecord = zeros(p,n);
% nldyn = zeros(p,n);
% pk = zeros(p,n);
% pkrecord = zeros(p,n);
EXKrec = zeros(p,n-1);
% exkrec = zeros(p,n-1);
% EYKrec = zeros(p-1,n-1);
% eykrec = zeros(p-1,n-1);
% sk = zeros(p-1,n-1);
% skrecord = zeros(p-1,n-1);

for test=1:NN
    
    xtrue = zeros(p,n);
    xtrue(:,1) = mvnrnd(x0,P0);
    
    for k=2:n
        W = mvnrnd(zeros(6,1),Q)';
        xtrue(:,k) = correct(NLdyn(xtrue(:,k-1),u,W));
    end
    
    xk = zeros(p,n);
    xk(:,1) = x0;
    Ppkm1 = diag([100 100 2 1000 1000 2]);
    
    for k=2:n
        % Lookup full state to linearize about
%         xig = xk(1,k-1); etag = xk(2,k-1);
        thetag = xk(3,k-1);
%         xia = xk(4,k-1); etaa = xk(5,k-1);
        thetaa = xk(6,k-1);
        
        % Jacobian estimate
        Fkm1 = Fpred(vg,thetag,va,thetaa);
        
        % Estimation-Error Covariance
        Pmk = Fkm1*Ppkm1*Fkm1' + Q;
        
        % State Estimate
        xmk = correct(NLdyn(xk(:,k-1),u,zeros(6,1)));
        % record NL state estimate
%         nldyn(:,k) = xmk;
        
        % Jacobian estimate
        Hk = Hpred(xmk(1),xmk(2),xmk(4),xmk(5));
        
        % Kalman Gain
        Kk  = Pmk*Hk' / (Hk*Pmk*Hk' + R);
        
        % Nonlinear Measurement Innovation
        ymkp1 = NLmeas(xmk);
        ymkp1(1) = wrapToPi(ymkp1(1)); ymkp1(3) = wrapToPi(ymkp1(3));
        
        eyk = [-angdiff(ydata(1,k),ymkp1(1)); ...
            ydata(2,k)-ymkp1(2); ...
            -angdiff(ydata(3,k),ymkp1(3)); ...
            ydata(4:5,k)-ymkp1(4:5)];
%         eykrec(:,k-1) = eyk;
        
        % NIS
%         epyk(k-1) = NIS(eyk,Hk,Pmk,R);
        
        % Update State Estimate
        xk(:,k) = correct(xmk + Kk*eyk);
        
        % Update Estimation-Error Covariance
        Ppkm1 = (I - Kk*Hk)*Pmk;
%         pk(:,k) = [2*sqrt(Ppk(1,1)) 2*sqrt(Ppk(2,2)) 2*sqrt(Ppk(3,3))...
%             2*sqrt(Ppk(4,4)) 2*sqrt(Ppk(5,5)) 2*sqrt(Ppk(6,6))]';
        
%         sk(:,k-1) = 2*sqrt(diag(Hk*Pmk*Hk' + R));
        
        
        
        % NEES
        diff = [xtrue(1:2,k)-xk(1:2,k); wrapToPi(xtrue(3,k)-xk(3,k)); ...
            xtrue(4:5,k)-xk(4:5,k); wrapToPi(xtrue(6,k)-xk(6,k))];
        nees(k-1) = nees(k-1) + NEES(diff,Ppkm1);
        EXKrec(:,k-1) = EXKrec(:,k-1) + diff;
%         Ppkm1 = Ppk;
    end
%     EYKrec = EYKrec + eykrec;
%     EXKrec = EXKrec + exkrec;
%     xkrecord = xkrecord + xk;
%     pkrecord = pkrecord + pk;
%     skrecord = skrecord + sk;
%     NLdynrecord = NLdynrecord + nldyn;
    
%     nis = nis+epyk;
%     nees = nees+exk;
end


if NN>1
    % EYKrec = EYKrec/NN;
    EXKrec = EXKrec/NN;
    % xkrecord = xkrecord/NN;
    % pkrecord = pkrecord/NN;
    % skrecord = skrecord/NN;
    % NLdynrecord = NLdynrecord/NN;
    % nis = nis/NN;
    nees = nees/NN;

    Observation = [EXKrec; nees];


    %% NEES
    alpha = 0.05;
    r1 = chi2inv(alpha/2,NN*6)/NN;
    r2 = chi2inv(1- alpha/2,NN*6)/NN;

    B = nees>r2;

    % Reward = r2 - nees.*B;
    % Reward = Reward + nees.*~B - r1^2;
    Reward = sum( (r2 - nees.*B) + (nees.*~B - r1^2) );
    %% NIS

    % r1 = chi2inv(alpha/2,NN*5)/NN;
    % r2 = chi2inv(1- alpha/2,NN*5)/NN;
    IsDone = 1;
end
% disp(Observation(:,1:10))
end