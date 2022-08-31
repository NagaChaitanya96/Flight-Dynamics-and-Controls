clc; close all; format shortG;
clear all;
global Thrust V c b W CD0 CDalp e CD_deltaE S CL0 CLalp CLq CLalp_dot CLdeltaE Cm0 Cmalp Cmq Cmalp_dot CmdeltaE
global Xu  Zu  U0 theta0 Mu Zw Xw Mw Zw_dot Mw_dot Zq Mq CD_0 CL_0 ef j
y04=zeros(1,4);
y02=zeros(1,2);
%% Aircraft Specifications
c=1.211 ;
b=10.47;
AR=8.8;
S=12.47;
m=750;
Ixx=873;
Iyy=907;
Izz=1680;
Ixz=1144;
Thrust=1136;
g=9.81;
W=m*g;
% Initial data for simulation
H=2000;
Vinf=60 ;
V1=60;       
V=V1;
%% Longitudinal data
CD0=0.036;
CDalp=0.041;
e=0.8;
CD_deltaE=0.026;
CL0=0.365;
CLalp=4.2;
CLq=17.3;
CLalp_dot=8.3;
CLdeltaE=0.26;
Cm0=0.05;
Cmalp=-.59;
Cmq=-9.3;
Cmalp_dot=-4.3;
CmdeltaE=-1.008;

%% Finding CL,CD,Cm values at trim condition
% Sea level conditions
rho0=1.225;
T0=288.15;
a=-6.5*10^-3;    % Lapse rate
T=T0+a.*H;       % Temperature at given altidude
% gradient layer
rho=rho0*(T/T0).^4.25588;
qinf=0.5*rho.*V.^2;
% At trim conditions L=W; T=D; Cm=0; theta=0 =>q=0
CL=(W/S)./(qinf);
CD=(Thrust/S)./(qinf);
Cm=0;
%% Calculation of stability derivatives
[Xu,Zu,U0,theta0,Mu,Zw,Xw,Mw,Zw_dot,Mw_dot,Zq,Mq,CD_0,CL_0]=longitudinal_derivatives(rho,AR,m,Iyy);
j=1;
Xderiv=[Xu;Xw]
Mderiv=[Mu;Mw;Mw_dot;Mq]
Zderiv=[Zu;Zw;Zw_dot;Zq]

Amatrix=[Xu(j)      Xw(j)        0      -g*cos(theta0)
    Zu(j)./(1-Zw_dot(j))         Zw(j)./(1-Zw_dot(j))       (U0(j)+Zq(j))./(1-Zw_dot(j))        -g*sin(theta0)./(1-Zw_dot(j))
    Mu(j)+((Mw_dot(j)*Zu(j))./(1-Zw_dot(j)))        Mw(j)+((Mw_dot(j)*Zw(j))./(1-Zw_dot(j)))        Mq(j)+((U0(j)+Zq(j))*Mw_dot(j))./(1-Zw_dot(j))      -g*sin(theta0)*Mw_dot(j)./(1-Zw_dot(j))
    0        0       1       0];
EigenA4=eig(Amatrix)


Ashort=[Zw(j)./(1-Zw_dot(j)) (U0(j)+Zq(j))./(1-Zw_dot(j))
    Mw(j)+(Mw_dot(j)*Zw(j))./(1-Zw_dot(j)) Mq(j)+((U0(j)+Zq(j))*Mw_dot(j))./(1-Zw_dot(j))];
EigShort=eig(Ashort)


Along=[Xu(j) -g*cos(theta0)
    -Zu(j)/(U0(j)+Zq(j)) g*sin(theta0)/(U0(j)+Zq(j))];
EigLong=eig(Along)
%% PLots of Eigen Values
figure
plot(EigenA4,'D','linewidth',1,'markeredgecolor','r','MarkerFacecolor','w')
xlabel('\zeta*\omega_n')
ylabel('Damping frequency')
hold on
plot(EigShort,'D','linewidth',1,'MarkerEdgecolor','b')
plot(EigLong,'D','MarkerEdgecolor','k')
legend('Fourth order system','Short Period mode','Long Period Mode')
title('Eigen Value plots for all Modes')
grid minor
axis equal

figure
plot(EigenA4,'D','linewidth',1,'markeredgecolor','r')
xlabel('\zeta*\omega_n')
ylabel('Damping frequency')
hold on
plot(EigShort,'D','MarkerEdgecolor','b','linewidth',1)
legend('Fourth order system','SHort Period Mode')
title('4rth order vs Short period ')
grid minor
axis auto

figure
plot(EigenA4,'D','linewidth',1,'markeredgecolor','r')
xlabel('\zeta*\omega_n')
ylabel('Damping frequency')
hold on
plot(EigLong,'D','MarkerEdgecolor','k','linewidth',1)
legend('Fourth order system','Long Period Mode')
title('4rth order vs long period mode')
grid minor
axis auto




