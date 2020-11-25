% ASEN 5044 Final Project
% Due: Friday, December 11, 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Housekeeping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

ImportCoop;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem 1

L = 0.5;
x_nom = [ 10;
           0;
          pi/2;
          -60;
           0;
         -pi/2];
u_nom = [   2;
         -pi/18;
           12;
          pi/25 ];
[ Atil_nom, Btil_nom, Gammatil_nom, Htil_nom ] = Problem_1(L,x_nom,u_nom);

%% Problem 2
dt=0.1;

Ftil = eye(6)+dt*Atil_nom;
Gtil = dt*Btil_nom;
Omegatil = dt*Gammatil_nom;
Htil=Htil_nom;

%% Problem 3

delx(:,1) = [0; 1; 0; 0; 0; 0.1];
delu = [0; 0; 0; 0];

Data=load('cooplocalization_finalproj_KFdata');
tvec=Data.tvec;
ydata=Data.ydata;

y_nom=ydata(:,2:end);

for k=1:length(tvec)
    
    x_nom(:,k+1)=x_nom(:,k)+delx(:,k);
    
    x(:,k+1) = x_nom(:,k+1)+Ftil*delx(:,k)+Gtil*delu;
    
    delx(:,k+1)=x(k+1)-x_nom(k+1);
    
    y(:,k) = y_nom(k)+Htil*delx(:,k+1);
    
    [ Atil_nom, Btil_nom, Gammatil_nom, Htil_nom ] = Problem_1(L,x_nom(:,k+1),u_nom);
    
    Ftil = eye(6)+dt*Atil_nom;
    Gtil = dt*Btil_nom;
    Omegatil = dt*Gammatil_nom;
    Htil = Htil_nom;
    
end


% Plotting states
figure;
subplot(6,1,1)
plot(tvec,x_nom(1,1:end-1),'k','LineWidth',3), hold on
xlabel('Time [s]','fontsize',22,'Interpreter','latex')
ylabel('$\xi_{g}$ [m]','fontsize',24,'Interpreter','latex')
set(gca, 'FontSize', 20)
grid minor;


subplot(6,1,2)
plot(tvec,x_nom(2,1:end-1),'k','LineWidth',3), hold on
xlabel('Time [s]','fontsize',22,'Interpreter','latex')
ylabel('$\eta_{g}$ [m]','fontsize',24,'Interpreter','latex')
set(gca, 'FontSize', 20)
grid minor;
     

subplot(6,1,3)
plot(tvec,x_nom(3,1:end-1),'k','LineWidth',3), hold on
xlabel('Time [s]','fontsize',22,'Interpreter','latex')
ylabel('$\theta_{g}$ [rad]','fontsize',24,'Interpreter','latex')
set(gca, 'FontSize', 20)
grid minor;


subplot(6,1,4)
plot(tvec,x_nom(4,1:end-1),'k','LineWidth',3), hold on
xlabel('Time [s]','fontsize',22,'Interpreter','latex')
ylabel('$\xi_{a}$ [m]','fontsize',24,'Interpreter','latex')
set(gca, 'FontSize', 20)
grid minor;


subplot(6,1,5)
plot(tvec,x_nom(5,1:end-1),'k','LineWidth',3), hold on
xlabel('Time [s]','fontsize',22,'Interpreter','latex')
ylabel('$\eta_{a}$ [m]','fontsize',24,'Interpreter','latex')
set(gca, 'FontSize', 20)
grid minor;
     

subplot(6,1,6)
plot(tvec,x_nom(6,1:end-1),'k','LineWidth',3), hold on
xlabel('Time [s]','fontsize',22,'Interpreter','latex')
ylabel('$\theta_{a}$ [rad]','fontsize',24,'Interpreter','latex')
set(gca, 'FontSize', 20)
grid minor;


suptitle('States vs Time Linearized Approximate Dynamics Simulation')


% Plotting perturbations
figure;
subplot(6,1,1)
plot(tvec,delx(1,1:end-1),'k','LineWidth',3), hold on
xlabel('Time [s]','fontsize',22,'Interpreter','latex')
ylabel('$\xi_{g}$ [m]','fontsize',24,'Interpreter','latex')
set(gca, 'FontSize', 20)
grid minor;


subplot(6,1,2)
plot(tvec,delx(2,1:end-1),'k','LineWidth',3), hold on
xlabel('Time [s]','fontsize',22,'Interpreter','latex')
ylabel('$\delta \eta_{g}$ [m]','fontsize',24,'Interpreter','latex')
set(gca, 'FontSize', 20)
grid minor;
     

subplot(6,1,3)
plot(tvec,delx(3,1:end-1),'k','LineWidth',3), hold on
xlabel('Time [s]','fontsize',22,'Interpreter','latex')
ylabel('$\delta \theta_{g}$ [rad]','fontsize',24,'Interpreter','latex')
set(gca, 'FontSize', 20)
grid minor;


subplot(6,1,4)
plot(tvec,delx(4,1:end-1),'k','LineWidth',3), hold on
xlabel('Time [s]','fontsize',22,'Interpreter','latex')
ylabel('$\delta \xi_{a}$ [m]','fontsize',24,'Interpreter','latex')
set(gca, 'FontSize', 20)
grid minor;


subplot(6,1,5)
plot(tvec,delx(5,1:end-1),'k','LineWidth',3), hold on
xlabel('Time [s]','fontsize',22,'Interpreter','latex')
ylabel('$\delta \eta_{a}$ [m]','fontsize',24,'Interpreter','latex')
set(gca, 'FontSize', 20)
grid minor;
     

subplot(6,1,6)
plot(tvec,delx(6,1:end-1),'k','LineWidth',3), hold on
xlabel('Time [s]','fontsize',22,'Interpreter','latex')
ylabel('$\delta \theta_{a}$ [rad]','fontsize',24,'Interpreter','latex')
set(gca, 'FontSize', 20)
grid minor;


suptitle('Linearized approx perturbations vs Time')


% Plotting Measurments
figure;
subplot(5,1,1)
plot(tvec,y(1,:),'k','LineWidth',3), hold on
xlabel('Time [s]','fontsize',22,'Interpreter','latex')
ylabel('$\gamma_{ag}$ [rad]','fontsize',24,'Interpreter','latex')
set(gca, 'FontSize', 20)
grid minor;


subplot(5,1,2)
plot(tvec,y(2,:),'k','LineWidth',3), hold on
xlabel('Time [s]','fontsize',22,'Interpreter','latex')
ylabel('$\rho_{ga}$ [m]','fontsize',24,'Interpreter','latex')
set(gca, 'FontSize', 20)
grid minor;
     

subplot(5,1,3)
plot(tvec,y(3,:),'k','LineWidth',3), hold on
xlabel('Time [s]','fontsize',22,'Interpreter','latex')
ylabel('$\gamma_{ga}$ [rad]','fontsize',24,'Interpreter','latex')
set(gca, 'FontSize', 20)
grid minor;


subplot(5,1,4)
plot(tvec,y(4,:),'k','LineWidth',3), hold on
xlabel('Time [s]','fontsize',22,'Interpreter','latex')
ylabel('$\xi_{a}$ [m]','fontsize',24,'Interpreter','latex')
set(gca, 'FontSize', 20)
grid minor;


subplot(5,1,5)
plot(tvec,y(5,:),'k','LineWidth',3), hold on
xlabel('Time [s]','fontsize',22,'Interpreter','latex')
ylabel('$\eta_{a}$ [m]','fontsize',24,'Interpreter','latex')
set(gca, 'FontSize', 20)
grid minor;


suptitle('Approximate Linearized Model Data Simulation')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem 4

% Part a

% Part b

% Part c

%% Problem 5

% Part a

% Part b

% Part c

%% Problem 6


%%%%%%%%%%%%%%%%%%%%%%%%%%% Advanced Questions %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% AQ13


%% AQ14

