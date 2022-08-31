clear all;close all; clc ;format shortg
global Thrust Vinf c b W CD0 CDalp e CD_deltaE S CL0 CLalp CLq CLalp_dot CL_deltaE Cm0 Cmalp Cmq Cmalp_dot Cm_deltaE
global CY0 CY_beta CYp CYr CY_deltaR CY_deltaA Cl0 Cl_beta Clp Clr Cl_deltaR Cl_deltaA Cn0 Cnalp Cnp Cnr Cn_deltaR Cn_deltaA
global CZ_deltaE Cn_beta

global theta0 A4lat B4lat dt phiold Kp Ki Kd
prompt=('Time for the disturbance to end in seconds  ');
dt=input(prompt)
controller =([' Type 1 for PI','\n Type 2 for PID', '\n Type 3 for PD','\n']);
ef=input(controller)
%% Controller type
if ef==1
Kp=-0.57832;
Ki=-0.34019;
Kd=0;
% Kp=-4.773; Kd=0; Ki=-1.573;
elseif ef==2
Kd=-0.38784;
Ki=-0.6666
Kp=-1.017
end
Gphic=tf(pid(Kp,Ki,Kd))

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
Cn0=0;Cn_beta=0.061;Cnp=-0.015;Cnr=-0.091;Cn_deltaR=-0.049;Cn_deltaA=0;
CZ_deltaE=-CL_deltaE;

% Aircraft transfer function for pitch angle
g=9.81;
W=m*9.81;
rho0=1.225;
T0=288.15;
a=-6.5*10^-3;
T=T0+a.*h;
rho=rho0*(T/T0).^4.25588;
Qinf=0.5.*rho*Vinf^2;

% Finding values
CL=(W/S)./Qinf;
CD1=(Thrust/W)./Qinf;
Cm=0;
CD2=CD0+(1/(pi*AR*e))*CL^2;

%% From stability equations finding trim angles in radians

A=[CLalp CLq CL_deltaE;CDalp 0 CD_deltaE   ;Cmalp Cmq Cm_deltaE];
dcl=CL -CL0;
B=[dcl;CD1-CD0;Cm-Cm0];
x=A\B;
alpha_t=x(1);   q_t=x(2);     deltaE_t=x(3);
% CD=CD0+CDalp.*abs(alpha)+CD_deltaE.*abs(alpha);
% [Xu,Xw,Zu,Zw,Zw_dot,Zq,Zalp,Zalp_dot,Mu,Mq,Mw,Malp,Mw_dot,Malp_dot,X_deltaE,Z_deltaE,M_deltaE,Cn_beta,Y_beta,Yp,Yr,N_beta,Np,Nr,L_beta,Lp,Lr,Y_deltaA,Y_deltaR,N_deltaA,N_deltaR,L_deltaA,L_deltaR,U0,theta0,CD_0,CL_0]=derivatives(rho,AR,m,Iyy,Ixx,Izz);
[Xu,Xw,Zu,Zw,Zw_dot,Zq,Zalp,Zalp_dot,Mu,Mq,Mw,Malp,Mw_dot,Malp_dot,X_deltaE,Z_deltaE,M_deltaE,Y_beta,Yp,Yr,N_beta,Np,Nr,L_beta,Lp,Lr,Y_deltaA,Y_deltaR,N_deltaA,N_deltaR,L_deltaA,L_deltaR,U0,theta0,CD_0,CL_0]=derivatives(rho,AR,m,Iyy,Ixx,Izz);

%% 4rth order LAteral system
Yv=Y_beta/U0;
Lv=L_beta/U0;
Nv=N_beta/U0;
c1=Ixz/Ixx;
c2=Ixz/Izz;
c3=c1*c2;
A4lat=[Yv Yp Yr-U0 g*cos(theta0)
    (Lv+Nv*c1)/(1-c3) (Lp+Np*c1)/(1-c3) (c1*Nr+Lr)/(1-c3)   0
    (Nv+Lv*c2)/(1-c3) (Np+Lp*c2)/(1-c3) (Lr*c2+Nr)/(1-c3)   0
    0 1 0 0];
B4lat=[0 Y_deltaR
    (L_deltaA+N_deltaA*c1)/(1-c3) (L_deltaR+N_deltaR*c1)/(1-c3)
    (N_deltaA+L_deltaA*c2)/(1-c3) (N_deltaR+L_deltaR*c2)/(1-c3)
    0 0];

B4aileron=B4lat(:,1);
C4lat=[0 0 0 1];
D4lat=0;

d_SSlat=ss(A4lat,B4aileron,C4lat,D4lat);
[num_phi,den_phi]=ss2tf(A4lat,B4aileron,C4lat,D4lat);
Gphi4=tf(num_phi,den_phi);
rlocus(Gphi4)

%% Simulation

t=0:01:60;
y0=[0 0 0 0.0];
deltaE_t=0;
phi_trim=0;
phiold=phi_trim;
[t,La]=ode45(@PID_roll_D,t,y0);
[t,LaN]=ode45(@dynamics_rol1,t,y0);
%% Plots
figure
for i=4
    plot(t,LaN(:,i),'linewidth',0.7)
    hold on
    plot(t,La(:,i))
end

xlabel('time (s)')
ylabel('phi (rad)')
if ef==1
title('Roll Dynamics with and without PI controller ')
legend('Actual system','With PI controller')
elseif ef==2
 title('Roll Dynamics with and without PID controller ')
legend('Actual system','With PID controller')   
end

figure
Gtotal=feedback(Gphi4*Gphic,1)
rlocus(Gtotal)
%% Functions for Dynamics of Aircraft
%% 1. Without Controller
function [dLaNdt]=dynamics_rol1(t,y)
global A4lat  dt
g=9.81;
dLaNdt=zeros(4,1);
v= y(1);
p=y(2) ;
r=y(3);
phi=y(4);
if  t>=0 && t<=dt
    phi=phi+0.5236;
end
dLaNdt(1)=A4lat(1,1)*v+ A4lat(1,2)*p + A4lat(1,3)*r+ A4lat(1,4)*phi; %+B4lat(1,1)*deltaA;
dLaNdt(2)=A4lat(2,1)*v+ A4lat(2,2)*p+ A4lat(2,3)*r +A4lat(2,4)*phi; %+B4lat(2,1)*deltaA;
dLaNdt(3)=A4lat(3,1)*v+A4lat(3,2)*p+ A4lat(3,3)*r+A4lat(3,4)*phi; %+B4lat(3,1)*deltaA;
dLaNdt(4)=A4lat(4,2)*p; %+B4lat(4,1)*deltaA;
end

%% 2. With Controller
function [dLadt]=PID_roll_D(t,y) 
global A4lat dt phiold B4lat Kp Ki Kd
dLadt=zeros(4,1);
% phi0=0; deltaAmax=0.6; deltaAmin=-0.6;
v= y(1); p=y(2) ; r=y(3); phi=y(4);
if  t>=1 && t<=dt
    phi=phi+0.5236;
end
phi_trim=0;
% if t>2
%     delta_aileron=Kp*(phi_trim-phi)-Kd*(p)+Ki*(phiold-phi);
% else
%     delta_aileron=0;
% end
delta_aileron=Kp*(phi_trim-phi)-Kd*(p)+Ki*(phiold-phi);
deltaA=delta_aileron;

dLadt(1)=A4lat(1,1)*v+ A4lat(1,2)*p + A4lat(1,3)*r+ A4lat(1,4)*phi +B4lat(1,1)*deltaA;
dLadt(2)=A4lat(2,1)*v+ A4lat(2,2)*p+ A4lat(2,3)*r +A4lat(2,4)*phi +B4lat(2,1)*deltaA;
dLadt(3)=A4lat(3,1)*v+A4lat(3,2)*p+ A4lat(3,3)*r+A4lat(3,4)*phi +B4lat(3,1)*deltaA;
dLadt(4)=A4lat(4,2)*p +B4lat(4,1)*deltaA;
phiold=phi;
end

