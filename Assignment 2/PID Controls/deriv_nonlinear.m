function [Xu,Xw,Zu,Zw,Zw_dot,Zq,Zalp,Zalp_dot,Mu,Mq,Mw,Malp,Mw_dot,Malp_dot,X_deltaE,Z_deltaE,M_deltaE,Cn_beta,Y_beta,Yp,Yr,N_beta,Np,Nr,L_beta,Lp,Lr,Y_deltaA,Y_deltaR,N_deltaA,N_deltaR,L_deltaA,L_deltaR,Yv,Lv,Nv]=deriv_nonlinear(CL,CD)
global   c b   CDalp e CD_deltaE S  CLalp CLq CLalp_dot CL_deltaE  Cmalp Cmq Cmalp_dot Cm_deltaE
global  CY_beta CYp CYr CY_deltaR CY_deltaA  Cl_beta Clp Clr Cl_deltaR Cl_deltaA   Cnp Cnr Cn_deltaR Cn_deltaA m Ixx Iyy Izz V
rho=1.0581;
g=9.81;
AR=12.775;
% V=55;
Q=0.5.*rho.*V.^2;
k=1/(pi*AR*e);
CDu=0;
CLu=0;
Cmu=0;
% Longitudinal Derivatives
Xu=         -(CDu+2*CD).*(Q.*S)./(m*V);
Xw=         (Q.*S./(m*V)).*(CL-CDalp);
  
Zu=         -(CLu+2*CL).*(Q.*S./(m*V));
Zw=         -(Q.*S./(m.*V)).*(CLalp+CD);
Zw_dot=         -(Q.*S*c./(2*m.*V.^2))*(CLalp_dot);
Zq=             -(Q.*S*c./(2*m.*V))*(CLq);
Zalp=V*Zw;
Zalp_dot=V*Zw_dot;

Mu=         (Q.*S*c./(Iyy.*V)).*(Cmu);
Mq=         (Q.*S*c^2./(2.*Iyy.*V))*(Cmq);
Mw=             (Q.*S*c./(V.*Iyy))*(Cmalp); 
Malp=V*Mw;
Mw_dot=         (Q.*S*c^2./(2.*Iyy.*V.^2))*(Cmalp_dot);
Malp_dot=V*Mw_dot;

X_deltaE=-CD_deltaE*Q*S./m;
Z_deltaE=-CL_deltaE*Q*S/m;
M_deltaE=Cm_deltaE*(Q*S*c)/Iyy;

Cn_beta=0;
%Lateral Derivatives
CY0=0;CY_beta=-0.531;CYp=0.2;CYr=0.633;CY_deltaR=0.15;CY_deltaA=0;
Y_beta=Q*S*b*CY_beta/m;
Yp=Q*S*b*CYp/(2*m*V);
Yr=Q*S*b*CYr/(2*m*V);

N_beta=Q*S*b*Cn_beta./Izz;
Np=Q*S*b^2*Cnp./(2*Izz*V);
Nr=Q*S*b^2*Cnr./(2*Izz*V);



L_beta=Q*S*b*Cl_beta./Ixx;
Lp=Q*S*b^2*Clp./(2*Ixx*V);
Lr=Q*S*b^2*Clr./(2*Ixx*V);


Y_deltaA=Q*S*CY_deltaA./m;
Y_deltaR=Q*S*CY_deltaR./m;

N_deltaA=Q*S*b*Cn_deltaA./Izz;
N_deltaR=Q*S*b*Cn_deltaR./Izz;

L_deltaA=Q*S*b*Cl_deltaA./Ixx;
L_deltaR=Q*S*b*Cl_deltaR./Ixx;

Yv=Y_beta/V;
Lv=L_beta/V;
Nv=N_beta/V;
end
