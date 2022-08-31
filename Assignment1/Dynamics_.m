clc; close all; format shortG;
clear all;
global Thrust V c b W CD0 CDalp e CD_deltaE S CL0 CLalp CLq CLalp_dot CLdeltaE Cm0 Cmalp Cmq Cmalp_dot CmdeltaE
global Xu  Zu  CL U0 theta0 Mu Zw Xw Mw Zw_dot Mw_dot Zq Mq CD_0 CL_0 ef j
prompt='Select the following modes 1. short period ;  2.longperiod ';
ef=input(prompt);
y04=zeros(1,4);
y02=zeros(1,2);
if ef==1
    tspan=0:0.01:30;
elseif ef==2
    tspan=0:0.01:400;
end
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
V1=[35 60 85 ];       % Variation of speed +-25m/s
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
% at trim condition deltaE=0; q=0;
abs_alpha=CL./CLalp;
%% Calculation of stability derivatives
[Xu,Zu,U0,theta0,Mu,Zw,Xw,Mw,Zw_dot,Mw_dot,Zq,Mq,CD_0,CL_0]=longitudinal_derivatives(rho,AR,m,Iyy);
Xderivatives=[Xu;Xw]
Mderiv=[Mu;Mw;Mw_dot;Mq]
Zderiv=[Zu;Zw;Zw_dot;Zq]

for j=1:length(V)
    Amatrix=[Xu(j)      Xw(j)        0      -g*cos(theta0)
        Zu(j)./(1-Zw_dot(j))         Zw(j)./(1-Zw_dot(j))       (U0(j)+Zq(j))./(1-Zw_dot(j))        -g*sin(theta0)./(1-Zw_dot(j))
        Mu(j)+((Mw_dot(j)*Zu(j))./(1-Zw_dot(j)))        Mw(j)+((Mw_dot(j)*Zw(j))./(1-Zw_dot(j)))        Mq(j)+((U0(j)+Zq(j))*Mw_dot(j))./(1-Zw_dot(j))      -g*sin(theta0)*Mw_dot(j)./(1-Zw_dot(j))
        0        0       1       0];
    EigenA4=eig(Amatrix);
    Fourth_O_eig_matrix(:,j)=EigenA4;
    %% Solving Ode for fourth order matix
    [t,y]=ode45(@characteristics,tspan,y04);
    ynames={'U','W','Q','theta'};
    %% Short period & Damping ratio and frequency
    if ef==1
        Ashort=[Zw(j)./(1-Zw_dot(j)) (U0(j)+Zq(j))./(1-Zw_dot(j))
            Mw(j)+(Mw_dot(j)*Zw(j))./(1-Zw_dot(j)) Mq(j)+((U0(j)+Zq(j))*Mw_dot(j))./(1-Zw_dot(j))];
        EigShort=eig(Ashort);
        EigShort_all(:,j)=EigShort;
        R_sp=abs(real(EigShort(1)));
        S_sp=imag(EigShort(1));
        val=R_sp./abs(EigShort);
        Zeta_Sp=val(1);
        val2=abs(EigShort).*(sqrt(1-Zeta_Sp.^2));
        Omega_sp=val2(1);
        Timeperiod_sp=(2*pi)/Omega_sp;
        T_half=log(2)/R_sp;
        disp(["                 Short Period mode Characteristics are            "])
        disp([' FOr velocity  ',num2str( V)])
        disp(['1.Damping ratio=',num2str(Zeta_Sp)])
        disp(['2.Frequency=',num2str(Omega_sp)])
        disp(['3.T_Halftime=', num2str(T_half)])
        disp(['4.TimePeriod=', num2str(Timeperiod_sp)])
        %% solving ODE for short period mode
        [t,sh]=ode45(@ShortPeriod,tspan,y02);
        %% PLOTTING SP VS 4 ORDER
        figure
        subplot(2,1,1)
        plot(t,y(:,2),'linewidth',1)
        hold on
        plot(t,sh(:,1),'--','linewidth',1)
        xlabel('time')
        ylabel('W (m/s)')
        hold off
        legend('Fourth Order system','Short Period Approximation','location','best')
        subplot(2,1,2)
        plot(t,y(:,3),'linewidth',1)
        hold on
        plot(t,sh(:,2),'--','linewidth',1)
        xlabel('time')
        ylabel('q (rad/s)')
        legend('Fourth Order system','Short Period Approximation','location','best')
        %%  Long period mode
    elseif ef==2
        Along=[Xu(j) -g*cos(theta0)
            -Zu(j)/(U0(j)+Zq(j)) g*sin(theta0)/(U0(j)+Zq(j))];
        EigLong=eig(Along);
        EigLong_all(:,j)=EigLong;
        R_ph=abs(real(EigLong(1)));
        S_ph=imag(EigLong(1));
        val3=R_ph./abs(EigLong);
        Zeta_ph=val3(1);
        val4=abs(EigLong).*(sqrt(1-Zeta_ph));
        Omega_ph=val4(1);
        Timeperiod_ph=(2*pi)/Omega_ph;
        T_half_l=log(2)/R_ph;
        disp(["                 Phugiod mode Characteristics are            "])
        disp(['1.Damping ratio=',num2str(Zeta_ph)])
        disp(['2.Frequency=',num2str(Omega_ph)])
        disp(['3.T_Halftime=', num2str(T_half_l)])
        disp(['4.TimePeriod=', num2str(Timeperiod_ph)])
        
        %% Solving the differential equations for long period mode
        [t,Lo]=ode45(@LongPeriod,tspan,y02);
        %% plotting LP VS 4 order
        figure
        subplot(2,1,1)
        plot(t,y(:,1),'linewidth',1)
        hold on
        plot(t,Lo(:,1),'--','linewidth',1)
        xlabel('time')
        ylabel('U (m/s)')
        legend('Fourth Order system','Phugoid Approximation')
        hold off
        subplot(2,1,2)
        plot(t,y(:,4),'linewidth',1)
        hold on
        plot(t,Lo(:,2),'--','linewidth',1)
        xlabel('time')
        ylabel('theta (rad)')
        hold off
        legend('Fourth Order system','Phugoid Approximation')
    end
end
