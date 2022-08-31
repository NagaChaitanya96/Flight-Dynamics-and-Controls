clear all;close all; clc ;format shortg
global Thrust Vinf c b W CD0 CDalp e CD_deltaE S CL0 CLalp CLq CLalp_dot CL_deltaE Cm0 Cmalp Cmq Cmalp_dot Cm_deltaE
global CY0 CY_beta CYp CYr CY_deltaR CY_deltaA Cl0 Cl_beta Clp Clr Cl_deltaR Cl_deltaA Cn0 Cnalp Cnp Cnr Cn_deltaR Cn_deltaA
global CZ_deltaE
global A4long B4long theta0  ef dt Kp Ki Kd W S V
global  thetaold m Ixx Iyy Izz

prompt=('Time for the disturbance to end in seconds  ');
dt=input(prompt);
controller =([' Type 1 for PID','\n Type 2 for PD', '\n Type 3 for PI','\n']);
ef=input(controller);



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
%% 4rth order Longitudinal system
A4long=[Xu      Xw        0      -g*cos(theta0)
    Zu./(1-Zw_dot)         Zw./(1-Zw_dot)       (U0+Zq)./(1-Zw_dot)        -g*sin(theta0)./(1-Zw_dot)
    Mu+((Mw_dot*Zu)./(1-Zw_dot))        Mw+((Mw_dot*Zw)./(1-Zw_dot))        Mq+((U0+Zq)*Mw_dot)./(1-Zw_dot)      -g*sin(theta0)*Mw_dot./(1-Zw_dot)
    0        0       1       0];
B4long=[X_deltaE;Z_deltaE;M_deltaE+Mw_dot*Z_deltaE;0];
C4long=[0 0 0 1];
D4long=0;
d_SS=ss(A4long,B4long,C4long,D4long);
[num_theta,den_theta]=ss2tf(A4long,B4long,C4long,D4long)
Gtheta4=tf(num_theta,den_theta)
damp(Gtheta4)
% figure
% rlocus(Gtheta4)

%% Simulation
t=0:0.001:400;
y0=[0 0 0 0];
deltaE_t=0;
theta_trim=0;
thetaold=theta_trim;
%% Controller type
if ef==1
    Kp=-3.829;
    Kd=-.6701;
    Ki=-5.47;
    [t,yn]=ode45(@dynamics_theta,t,y0);
    [t,y]=ode45(@PID_theD,t,y0);
    
    GthetaC=tf(pid(Kp,Ki,Kd))
    figure;clf
    plot(t,yn(:,4),'linewidth',1.1)
    hold on
    plot(t,y(:,4),'linewidth',1.1)
    xlabel('time (s)')
    ylabel('theta (rad)')
    legend('theta w/o controller','theta PID controller')
    subtitle('Effect of PID on Theta for nonlinear modal')
    
elseif ef==2
    Kp=-0.060;
    Kd=-0.518;
    Ki=0;
    [t,yn]=ode45(@dynamics_theta,t,y0);
    [t,y]=ode45(@PID_theD,t,y0);
    GthetaC=tf(pid(Kp,Ki,Kd))
    figure;clf
    
    plot(t,yn(:,4),'linewidth',1.1)
    hold on
    plot(t,y(:,4),'linewidth',1.1)
    xlabel('time (s)')
    ylabel('theta (rad)')
    legend('theta w/o controller','theta PD controller')
    subtitle('Effect of PD on Theta for nonlinear modal')
    
elseif ef==3
    Ki=-0.02862983;
    Kp=Ki*15;
    Kd=0;
% Kp=-0.21; Ki=-0.08275; Kd=0;
    [t,yn]=ode45(@dynamics_theta,t,y0);
    [t,y]=ode45(@PID_theD,t,y0);
    GthetaC=tf(pid(Kp,Ki,Kd))
    figure;clf
    plot(t,yn(:,4),'linewidth',1.1)
    hold on
    plot(t,y(:,4),'linewidth',1.1)
    xlabel('time (s)')
    ylabel('theta (rad)')
    legend('theta w/o controller','theta PI controller')
    subtitle('Effect of PI on Theta for nonlinear modal')
    
end
figure
Gtotal=feedback(Gtheta4*GthetaC,1)
rlocus(Gtotal)
%%                          Functions for Dynamics of Aircraft

%% 1. Without Controller
function [dyndt]=dynamics_theta(t,yn) %deltaE_t
global A4long  dt
g=9.81;
dyndt=zeros(4,1);
u= yn(1); w=yn(2) ; q=yn(3); theta=yn(4);
if  t>=1 && t<=dt
    theta=theta+0.087266;
end
dyndt(1)=A4long(1,1)*u+ A4long(1,2)*w + A4long(1,3)*q+ A4long(1,4)*theta; %+B4long(1)*deltaE;
dyndt(2)=A4long(2,1)*u+ A4long(2,2)*w+ A4long(2,3)*q +A4long(2,4)*theta; %+B4long(2)*deltaE;
dyndt(3)=A4long(3,1)*u+A4long(3,2)*w+ A4long(3,3)*q+A4long(3,4)*theta ;%+B4long(3)*deltaE;
dyndt(4)=A4long(4,3)*q; %+ B4long(4)*deltaE;

end


%% 2. With Controller

function [dydt]=PID_theD(t,y) %deltaE_t
global  dt thetaold   Kp Ki Kd W S c b  CLalp CL_deltaE Cmalp Cm_deltaE CL0 Cm0 CLq CD0 CDalp CD_deltaE Cmq V
%m Ixx Iyy Izz
dydt=zeros(4,1);
u= y(1); w=y(2) ; q=y(3); theta=y(4);
if  t>=1 && t<=dt
    theta=theta+0.087266;
end
U=55;
h=1500;
U=U+u;
V=norm(U,w);
rho=1.0581;
CL_trim=2*(W/S)/(rho*V^2);
A=[CLalp CL_deltaE;Cmalp Cm_deltaE];
B=[CL_trim-CL0 ; -Cm0];

angles_trim=A\B;

alpha_trim=angles_trim(1);
deltaE_trim=angles_trim(2);
theta_trim=0;

% if t>5
%     delta_elevator=Kp*(theta_trim-theta)-Kd*(q)+Ki*(thetaold-theta);
%     deltaE=deltaE_trim+delta_elevator;
% else
%     deltaE=deltaE_trim;
% end
delta_elevator=Kp*(theta_trim-theta)-Kd*(q)+Ki*(thetaold-theta);
deltaE=deltaE_trim+delta_elevator;
alpha=alpha_trim+(w/V);

CL=CL0+CLalp*(alpha)+CLq*q*c/(2*V)+CL_deltaE*(deltaE);
CD=CD0+CDalp*abs(alpha)+CD_deltaE*abs(deltaE);
Cm=Cm0+Cmalp*(alpha)+Cmq*(q*c/(2*V))+Cm_deltaE*deltaE;
CD2=CD0+(c/b)*CL^2;

[Xu,Xw,Zu,Zw,Zw_dot,Zq,Zalp,Zalp_dot,Mu,Mq,Mw,Malp,Mw_dot,Malp_dot,X_deltaE,Z_deltaE,M_deltaE,Cn_beta,Y_beta,Yp,Yr,N_beta,Np,Nr,L_beta,Lp,Lr,Y_deltaA,Y_deltaR,N_deltaA,N_deltaR,L_deltaA,L_deltaR]=deriv_nonlinear(CL,CD);
g=9.81;
U0=V;
Anl=[Xu      Xw        0      -g*cos(theta)
    Zu./(1-Zw_dot)         Zw./(1-Zw_dot)       (U0+Zq)./(1-Zw_dot)        -g*sin(theta)./(1-Zw_dot)
    Mu+((Mw_dot*Zu)./(1-Zw_dot))        Mw+((Mw_dot*Zw)./(1-Zw_dot))        Mq+((U0+Zq)*Mw_dot)./(1-Zw_dot)      -g*sin(theta)*Mw_dot./(1-Zw_dot)
    0        0       1       0];
Bnl=[X_deltaE;Z_deltaE;M_deltaE+Mw_dot*Z_deltaE;0];


dydt(1)=Anl(1,1)*u+ Anl(1,2)*w + Anl(1,3)*q+ Anl(1,4)*theta    +Bnl(1)*deltaE;
dydt(2)=Anl(2,1)*u+ Anl(2,2)*w+ Anl(2,3)*q +Anl(2,4)*theta     +Bnl(2)*deltaE;
dydt(3)=Anl(3,1)*u+Anl(3,2)*w+ Anl(3,3)*q+Anl(3,4)*theta   +Bnl(3)*deltaE;
dydt(4)=Anl(4,3)*q       + Bnl(4)*deltaE;
thetaold=theta;
end