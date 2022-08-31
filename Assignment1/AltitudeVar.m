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
H=[2000 3500 5000];
Vinf=60 ;
V=60;
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
%% Finding alpha, deltaE, q from the CL,CD,Cm, and from Longitudinal model equations
for i=1:length(H)
    A=[CLalp CLdeltaE;Cmalp CmdeltaE];
    dcl(i)=CL(i) -CL0  ;   % delta cl
    B=[dcl(i);Cm-Cm0];
    x=A\B;
    X(:,i)=x;
    alpha=X(1,:); q=0;
    deltaE=X(2,:);
    CD1=CD0+CDalp.*abs(alpha)+CD_deltaE.*abs(alpha);
end
%% Calculation of stability derivatives
[Xu,Zu,U0,theta0,Mu,Zw,Xw,Mw,Zw_dot,Mw_dot,Zq,Mq,CD_0,CL_0]=longitudinal_derivatives(rho,AR,m,Iyy);
Xderivatives=[Xu;Xw]
Mderiv=[Mu;Mw;Mw_dot;Mq]
Zderiv=[Zu;Zw;Zw_dot;Zq]

for j=1:length(H)
    Amatrix=[Xu(j)      Xw(j)        0      -g*cos(theta0)
        Zu(j)./(1-Zw_dot(j))         Zw(j)./(1-Zw_dot(j))       (U0+Zq(j))./(1-Zw_dot(j))        -g*sin(theta0)./(1-Zw_dot(j))
        Mu(j)+((Mw_dot(j)*Zu(j))./(1-Zw_dot(j)))        Mw(j)+((Mw_dot(j)*Zw(j))./(1-Zw_dot(j)))        Mq(j)+((U0+Zq(j))*Mw_dot(j))./(1-Zw_dot(j))      -g*sin(theta0)*Mw_dot(j)./(1-Zw_dot(j))
        0        0       1       0];
    EigenA4=eig(Amatrix);
    Fourth_O_eig_matrix(:,j)=EigenA4;
    %% Short period & Damping ratio and frequency
    Ashort=[Zw(j)./(1-Zw_dot(j)) (U0+Zq(j))./(1-Zw_dot(j))
        Mw(j)+(Mw_dot(j)*Zw(j))./(1-Zw_dot(j)) Mq(j)+((U0+Zq(j))*Mw_dot(j))./(1-Zw_dot(j))];
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
    
    Along=[Xu(j) -g*cos(theta0)
        -Zu(j)/(U0+Zq(j)) g*sin(theta0)/(U0+Zq(j))];
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
end
Xderivatives=[Xu;Xw]
Mderiv=[Mu;Mw;Mw_dot;Mq]
Zderiv=[Zu;Zw;Zw_dot;Zq]
%% PLotting Locus of Eigen values
% Variation in Altitude
figure
axis equal
grid minor
for j=1:length(H)
    hold on
    plot(Fourth_O_eig_matrix(:,j),'d','linewidth',1)
    title('Plotting all Eigen values locus for variation in Height')
    xlabel( '{\zeta\omega_n}')
    ylabel('{\omega_d}')
end

legend('H=2000m','H=3500m','H=5000 m','location','best')

figure
axis equal
grid minor
for j=1:length(H)
    hold on
    plot(EigShort_all(:,j),'d','linewidth',1)
    title('Eigen Values for Short period mode')
    xlabel( '{\zeta\omega_n}')
    ylabel('{\omega_d}')
    
end

legend('H=2000m','H=3500m','H=5000 m','location','best')

figure
axis auto
grid minor
for j=1:length(H)
    hold on
    plot(EigLong_all(:,j),'d','linewidth',1)
    title('Eigen Values for Phugoid mode')
    xlabel( '{\zeta\omega_n}')
    ylabel('{\omega_d}')
end

legend('H=2000m','H=3500m','H=5000 m','location','best')
%% Limitation
WL=W/S;
TW=T./W;
AE=CL./CD;
Fourth_O_eig_matrix(:,3)
EigShort_all(:,3)
EigLong_all(:,3)
figure
plot(Fourth_O_eig_matrix(:,3),'d')
hold on
plot(EigShort_all(:,3),'d')
plot(EigLong_all(:,3),'d')
legend('4order','sp','lp')
grid minor
% axis equal
%% Characteristics of fourth order vs reduced order model
gta= abs(real(Fourth_O_eig_matrix(:,3)));
eta=imag(Fourth_O_eig_matrix(:,3));
ab=(gta.^2+eta.^2).^0.5;
val_all=gta./ab;
% Zeta
g1=val_all(1);
g2=val_all(3);
% Omega_N
val_4sp=ab(1).*(sqrt(1-g1.^2));
val_4long=ab(3).*(sqrt(1-g2.^2));

Omega_4sp=val_4sp;
Timeperiod_4sp=(2*pi)/Omega_4sp;
T_half_4sp=log(2)/gta(1);


Omega_4long=val_4long;
Timeperiod_4long=(2*pi)/Omega_4long;
T_half_4long=log(2)/gta(3);

T4sp=[g1;Timeperiod_4sp;Omega_4sp;T_half_4sp];
Tsp=[Zeta_Sp;Timeperiod_sp;Omega_sp;T_half];
variables={'zeta';'Timeperiod';'Omega_damping';'T Half'};
T4lp=[g2;Timeperiod_4long;Omega_4long;T_half_4long];
Tlp=[Zeta_ph;Timeperiod_ph;Omega_ph;T_half_l];

error_sp=abs([T4sp-Tsp]);
error_lp=abs([T4lp-Tlp]);
T1= table(variables,T4sp,Tsp,error_sp);
T2=table(variables,T4lp,Tlp,error_lp);
disp([])
disp(['Wingloading= ', num2str(WL)])
disp(['ThrusttoWeight_ratio= ', num2str(TW)])
disp(['Aerodynamics_efficiency= ', num2str(AE)])
disp([T1])
disp([T2])
