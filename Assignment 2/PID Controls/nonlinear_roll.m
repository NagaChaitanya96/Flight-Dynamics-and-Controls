clear all;close all; clc ;format shortg
global Thrust Vinf c b W CD0 CDalp e CD_deltaE S CL0 CLalp CLq CLalp_dot CL_deltaE Cm0 Cmalp Cmq Cmalp_dot Cm_deltaE
global CY0 CY_beta CYp CYr CY_deltaR CY_deltaA Cl0 Cl_beta Clp Clr Cl_deltaR Cl_deltaA Cn0 Cnalp Cnp Cnr Cn_deltaR Cn_deltaA
global CZ_deltaE Cn_beta Ixz Ixx Izz
global theta0 A4lat B4lat dt phiold Kp Ki Kd
prompt=('Time for the disturbance to end in seconds  ');
dt=input(prompt);
controller =([' Type 1 for PI','\n Type 2 for PID', '\n ']);
ef=input(controller);

% s=tf('s')
%  -0.062576* (s+52.41) *(s+2.576)/s
%% Controller type
if ef==1
Kp=-0.57832;
Ki=-0.34019;
Kd=0;
% Kp=-4.773; Kd=0; Ki=-1.573;
elseif ef==2
% Kd=-0.38784;
% Ki=-0.6666
% Kp=-1.017
Kp=-3.441; Kd=-0.06258; Ki=-8.448;
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
title('Roll Dynamics with and without PI controller for non linear model')
legend('Actual system','With PI controller')
elseif ef==2
 title('Roll Dynamics with and without PID controller for non linear model ')
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
global  dt phiold   Kp Ki Kd W S c b tA4lat B4la m
global Blat_nl   Alat_nl CLalp CL_deltaE Cmalp Cm_deltaE CL0 Cm0 CLq CD0 CDalp CD_deltaE Cmq V
global CY0 CY_beta CYp CYr CY_deltaR CY_deltaA Cl0 Cl_beta Clp Clr Cl_deltaR Cl_deltaA Cn0 Cnalp Cnp Cnr Cn_deltaR Cn_deltaA Cn_beta

dLadt=zeros(4,1);
v= y(1); p=y(2) ; r=y(3); phi=y(4);
if  t>=1 && t<=dt
    phi=phi+0.5236;
end
U=55;
h=1500;
V=sqrt(U^2+v^2);
rho=1.0581;
CL_trim=2*(W/S)/(rho*V^2);
CY_trim=2*(W*sin(phi)/S)/(1.225*V^2);

beta_trim=0;
deltaA_trim=0;

phi_trim=0;
% % if t>5
%     delta_elevator=Kp*(theta_trim-theta)-Kd*(q)+Ki*(thetaold-theta);
%     deltaE=deltaE_trim+delta_elevator;
% else
%     deltaE=deltaE_trim;
% end
delta_aileron=Kp*(phi_trim-phi)-Kd*(p)+Ki*(phiold-phi);
% delta_aileron=Kp*(phi_trim-phi)-Kd*(p+r)+Ki*(phiold-phi);

deltaA=deltaA_trim+delta_aileron;
deltaR=0;
beta=beta_trim+v/U;

CY=CY0+CY_beta*beta+CYp*(p*b/(2*V))+CYr+(r*b/(2*V))+Cl_deltaR*(deltaR);
Cl=Cl0+Cl_beta*(beta)+Clp*(p*b/(2*V))+Clr*(r*b/(2*V))+Cl_deltaA*deltaA+Cl_deltaR*deltaR;
Cn=Cn0+Cn_beta*beta+Cnp*(p*b/(2*V))+Cnr*(r*b/(2*V))+Cn_deltaR*deltaR;

CL=CL_trim;
CD=CD0+(S/b^2)*CL^2;
% [Xu,Xw,Zu,Zw,Zw_dot,Zq,Zalp,Zalp_dot,Mu,Mq,Mw,Malp,Mw_dot,Malp_dot,X_deltaE,Z_deltaE,M_deltaE,Cn_beta,Y_beta,Yp,Yr,N_beta,Np,Nr,L_beta,Lp,Lr,Y_deltaA,Y_deltaR,N_deltaA,N_deltaR,L_deltaA,L_deltaR,Yv,Lv,Nv]=deriv_nonlinear(CL,CD);

[Y_beta,Yp,Yr,N_beta,Np,Nr,L_beta,Lp,Lr,Y_deltaA,Y_deltaR,N_deltaA,N_deltaR,L_deltaA,L_deltaR,Yv,Lv,Nv,k,c1,c2,c3]=nonlinear_deriv_lateral();

g=9.81;
U0=V;
theta0=0;
%% Lateral derivatives
rho=1.0581;
Cy0=0;Cy_beta=-0.531;Cyp=0.2;Cyr=0.633;Cy_deltaR=0.15;Cy_deltaA=0;
Q=0.5*rho*V^2;
Ybeta=Q*12.47*(-0.531)/700;
% Yv=Q*S*Cy_beta/(m*V);

Yp=0;
Yr=Q*12.47*15.47*0.633/(2*700*V);
Y_deltaA=0;
Y_deltaR=Q*12.47*0.15/700;
Yv=Ybeta/V;
Nv=N_beta/V;
Alat_nl=[Yv Yp Yr-U0 0;
    (Lv+Nv*c1)/(1-c3) (Lp+Np*c1)/(1-c3) (c1*Nr+Lr)/(1-c3)   0;
    (Nv+Lv*c2)/(1-c3) (Np+Lp*c2)/(1-c3) (Lr*c2+Nr)/(1-c3)   0;
    0 1 0 0];
Blat_nl=[0 Y_deltaR
    (L_deltaA+N_deltaA*c1)/(1-c3) (L_deltaR+N_deltaR*c1)/(1-c3)
    (N_deltaA+L_deltaA*c2)/(1-c3) (N_deltaR+L_deltaR*c2)/(1-c3)
    0 0];

dLadt(1)=Alat_nl(1,1)*v+ Alat_nl(1,2)*p + Alat_nl(1,3)*r+ Alat_nl(1,4)*phi    +Blat_nl(1,1)*deltaA;
dLadt(2)=Alat_nl(2,1)*v+ Alat_nl(2,2)*p+ Alat_nl(2,3)*r +Alat_nl(2,4)*phi     +Blat_nl(2,1)*deltaA;
dLadt(3)=Alat_nl(3,1)*v+Alat_nl(3,2)*p+ Alat_nl(3,3)*r+Alat_nl(3,4)*phi   +Blat_nl(3,1)*deltaA;
dLadt(4)=Alat_nl(4,2)*p;  + Blat_nl(4,1)*deltaA;
phiold=phi;

end

% lateral includes sideslip => v varies neglecting w,u variation.
% V=norm(U,v);


% A=[CY_beta CY_deltaA;Cn_beta Cn_deltaA]
% B=[CY_trim-CY0 ; -Cn0]

% angles_trim=A\B;
% 
% beta_trim=angles_trim(1)
% deltaA_trim=angles_trim(2)

% % 
% q=0;
% alpha=alpha_trim;
% deltaE=deltaE_trim;


% CL=CL0+CLalp*(alpha)+CLq*q*c/(2*V)+CL_d6eltaE*(deltaE);
% CD=CD0+CDalp*abs(alpha)+CD_deltaE*abs(deltaE);
% Cm=Cm0+Cmalp*(alpha)+Cmq*(q*c/(2*V))+Cm_deltaE*deltaE;
% CD2=CD0+(c/b)*CL^2;

% From symmetry of aircraft rudder at trim=0;