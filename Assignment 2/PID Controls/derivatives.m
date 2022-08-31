function [Xu,Xw,Zu,Zw,Zw_dot,Zq,Zalp,Zalp_dot,Mu,Mq,Mw,Malp,Mw_dot,Malp_dot,X_deltaE,Z_deltaE,M_deltaE,Y_beta,Yp,Yr,N_beta,Np,Nr,L_beta,Lp,Lr,Y_deltaA,Y_deltaR,N_deltaA,N_deltaR,L_deltaA,L_deltaR,U0,theta0,CD_0,CL_0]=derivatives(rho,AR,m,Iyy,Ixx,Izz)


global Thrust Vinf c b W CD0 CDalp e CD_deltaE S CL0 CLalp CLq CLalp_dot CL_deltaE Cm0 Cmalp Cmq Cmalp_dot Cm_deltaE
global CY0 CY_beta CYp CYr CY_deltaR CY_deltaA Cl0 Cl_beta Clp Clr Cl_deltaR Cl_deltaA Cn0 Cnalp Cnp Cnr Cn_deltaR Cn_deltaA
global  Cn_beta

U0=Vinf;
Q0=0.5.*rho.*U0.^2;
Q=Q0;
k=1/(pi*AR*e);
theta0=0;         %q*t
CL_0=W.*cos(theta0)./(Q0*S); % CL at alpha zero vec
CD_0=CD0+k*CL_0.^2; %vec
CDu=0;
CLu=0;
Cmu=0;
% Longitudinal Derivatives
Xu=         -(CDu+2*CD_0).*(Q0.*S)./(m*U0);
Xw=         (Q0.*S./(m*U0)).*(CL_0-CDalp);
  
Zu=         -(CLu+2*CL_0).*(Q0.*S./(m*U0));
Zw=         -(Q0.*S./(m.*U0)).*(CLalp+CD_0);
Zw_dot=         -(Q0.*S*c./(2*m.*U0.^2))*(CLalp_dot);
Zq=             -(Q0.*S*c./(2*m.*U0))*(CLq);
Zalp=U0*Zw;
Zalp_dot=U0*Zw_dot;

Mu=         (Q0.*S*c./(Iyy.*U0)).*(Cmu);
Mq=         (Q0.*S*c^2./(2.*Iyy.*U0))*(Cmq);
Mw=             (Q0.*S*c./(U0.*Iyy))*(Cmalp); 
Malp=U0*Mw;
Mw_dot=         (Q0.*S*c^2./(2.*Iyy.*U0.^2))*(Cmalp_dot);
Malp_dot=U0*Mw_dot;

X_deltaE=-CD_deltaE*Q0*S/m;
Z_deltaE=-CL_deltaE*Q0*S/m;
M_deltaE=Cm_deltaE*(Q0*S*c)/Iyy;

%Lateral Derivatives
Y_beta=Q0*S*CY_beta./m;
Yp=Q*S*b*CYp/(2*m*U0);
Yr=Q*S*b*CYr/(2*m*U0);

N_beta=Q0*S*b*Cn_beta./Izz;
Np=Q0*S*b^2*Cnp./(2*Izz*U0);
Nr=Q0*S*b^2*Cnr./(2*Izz*U0);



L_beta=Q0*S*b*Cl_beta./Ixx;
Lp=Q0*S*b^2*Clp./(2*Ixx*U0);
Lr=Q0*S*b^2*Clr./(2*Ixx*U0);


Y_deltaA=Q0*S*CY_deltaA./m;
Y_deltaR=Q0*S*CY_deltaR./m;

N_deltaA=Q0*S*b*Cn_deltaA./Izz;
N_deltaR=Q0*S*b*Cn_deltaR./Izz;

L_deltaA=Q0*S*b*Cl_deltaA./Ixx;
L_deltaR=Q0*S*b*Cl_deltaR./Ixx;
end
