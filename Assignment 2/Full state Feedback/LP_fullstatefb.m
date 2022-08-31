
clear all;close all; clc ;format shortg
global Thrust Vinf c b W CD0 CDalp e CD_deltaE S CL0 CLalp CLq CLalp_dot CL_deltaE Cm0 Cmalp Cmq Cmalp_dot Cm_deltaE
global CY0 CY_beta CYp CYr CY_deltaR CY_deltaA Cl0 Cl_beta Clp Clr Cl_deltaR Cl_deltaA Cn0 Cnalp Cnp Cnr Cn_deltaR Cn_deltaA
global CZ_deltaE error1
global A5long B5long theta0  ef dt Kp Ki Kd W S V Iter idx
global  thetaold m Ixx Iyy Izz T rho T0 rho0  AR  m Iyy Ixx Izz error alpha_trim q_trim deltaE_trim
global Alp_w_control a11 a12 a21 a22 ef nc11 nc12 nc21 nc22 ef

Kp=- 0.00071804;
Kd=-.0036445;
Ki=-3.5367e-05;
prompt=('Time for the disturbance to end in seconds  ');
dt=input(prompt);
%% GIven FLight Data
c=1.211;b=15.47;AR=19.9;S=12.47;m=700;
Ixx=1073;Iyy=907;Izz=1680;Ixz=1144;Thrust=800;
Vinf=55; h=1500;

%% Aerodynamic Derivatives
% Longitudinal
CD0=0.036;CDalp=0.061;e=0.9;CD_deltaE=0.026;
CL0=0.365;CLalp=5.2;CLq=17.3;CL_deltaE=0.5;CLalp_dot=0;
Cm0=0.05;Cmalp=-0.529;Cmq=-5.8;Cm_deltaE=-1.28;Cmalp_dot=0;

% Lateral-directional
CY0=0;CY_beta=-0.531;CYp=0.2;CYr=0.633;CY_deltaR=0.15;CY_deltaA=0;
Cl0=0;Cl_beta=-0.031;Clp=-0.3251;Clr=0.1836;Cl_deltaR=0.005;Cl_deltaA=-0.153;
Cn0=0;Cnalp=0.061;Cnp=-0.015;Cnr=-0.091;Cn_deltaR=-0.049;Cn_deltaA=0;
CZ_deltaE=-CL_deltaE;

% Aircraft transfer function for pitch angle
g=9.81;
W=m*9.81;
rho0=1.225;
T0=288.15;
a=-6.5*10^-3;
T=T0+a.*h;
rho=rho0*(T/T0).^4.25588
Qinf=0.5.*rho*Vinf^2
% Finding values
CL=(W/S)./Qinf
CD1=(Thrust/S)./Qinf %7.2797*1e-5 this is impossible
Cm=0;
CD2=CD0+(1/(pi*AR*e))*CL^2
CD = CD1+CD2

%% Trim results of angles , CL ,CD; CM =0
CL_trim=2*(W/S)/(rho*Vinf^2) % this formula can be applied at trim in steady level flight
CD_trim = 2*(Thrust/S)./(rho*Vinf^2) % same as known
A=[CLalp CLq CL_deltaE;CDalp 0 CD_deltaE   ;Cmalp Cmq Cm_deltaE];
dcl=CL_trim -CL0;
B=[dcl;CD_trim-CD0;0-Cm0];
x=A\B
%%
alpha=x(1);   q=x(2);     deltaE=x(3);
[alpha q deltaE]
CD=CD0+CDalp.*(alpha)+CD_deltaE.*abs(alpha);
[Xu,Xw,Zu,Zw,Zw_dot,Zq,Zalp,Zalp_dot,Mu,...
    Mq,Mw,Malp,Mw_dot,Malp_dot,X_deltaE,Z_deltaE,M_deltaE,...
    Y_beta,Yp,Yr,N_beta,Np,Nr,L_beta,Lp,Lr,Y_deltaA,....
    Y_deltaR,N_deltaA,N_deltaR,L_deltaA,L_deltaR,U0,...
    theta0,CD_0,CL_0]=derivatives(rho,AR,m,Iyy,Ixx,Izz)
%% Full state feedback By pole placement method
% taking U=-K*X and reference at trim r=0
%% Parameters taken for controller are 
% damping constant =0.3 
% Natural frequency =1.5 rad/s
% disturbance given for 20 seconds
% pitch autopilot
syms k1 k2
eqn1=-0.74123*k1-0.25917*k2==0.860438;
eqn2=-0.0156576*k2-2.5423*k1==2.31356;
Value=solve([eqn1,eqn2],[k1 k2])
% value=sym2cell(Value)
k1=double(Value.k1);
k2=double(Value.k2);

j=1;ef=2;
if ef==1
    tspan=0:0.01:30;
elseif ef==2
    tspan=0:0.01:400;
    y02=[0 0];
end
Alp_wo_control=[Xu(j) -g*cos(theta0)
    -Zu(j)/(U0(j)+Zq(j)) g*sin(theta0)/(U0(j)+Zq(j))]

nc11=Alp_wo_control(1,1);
nc12=Alp_wo_control(1,2);
nc21=Alp_wo_control(2,1);
nc22=Alp_wo_control(2,2);
[t,Lo_nc]=ode45(@Lp_without_c,tspan,y02);

Alp_w_control=[Xu-k1*X_deltaE -(g+k2*X_deltaE); (-Zu/U0+k1*Zu/U0) k2*Zu/U0]
a11=Alp_w_control(1,1);
a12=Alp_w_control(1,2);
a21=Alp_w_control(2,1);
a22=Alp_w_control(2,2);
[t,Lp]=ode45(@LP_fsfb,tspan,y02);
figure
% plot()
title('Comparision between FSFB and without controller for variation in theta=5 degrees')
subplot(2,1,1)
plot(t,Lo_nc(:,1),'linewidth',1)
hold on
plot(t,Lp(:,1),'linewidth',1)
xlabel('time')
ylabel('U (m/s)')
legend(' Phugoid without controller ','Phugoid with full state feedback controller')
title('Comparision between FSFB and without controller for variation in theta=5 degrees')

subplot(2,1,2)
plot(t,Lo_nc(:,2),'linewidth',1)
hold on
plot(t,Lp(:,2),'linewidth',1)
xlabel('time')
ylabel('theta (rad)')
legend(' Phugoid without controller ','Phugoid with full state feedback controller')
title('Comparision between FSFB and without controller for variation in theta=5 degrees')
