function [Xu,Zu,U0,theta0,Mu,Zw,Xw,Mw,Zw_dot,Mw_dot,Zq,Mq,CD_0,CL_0]=longitudinal_derivatives(rho,AR,m,Iyy)
global Thrust V c b W CD0 CDalp e CD_deltaE S CL0 CLalp CLq CLalp_dot CLdeltaE Cm0 Cmalp Cmq Cmalp_dot CmdeltaE
%.
U0=V;
Q0=0.5.*rho.*U0.^2;
k=1/(pi*AR*e);
theta0=0;         %q*t
CL_0=W.*cos(theta0)./(Q0*S); % CL at alpha zero vec
CD_0=CD0+k*CL_0.^2; %vec
CDu=0;
CLu=0;
Cmu=0;
% all are vectors
Xu=         -(CDu+2*CD_0).*(Q0.*S)./(m*U0);  
    
Zu=         -(CLu+2*CL_0).*(Q0.*S./(m*U0)); 

Mu=         (Q0.*S*c./(Iyy.*U0)).*(Cmu);

Zw=         -(Q0.*S./(m.*U0)).*(CLalp+CD_0);

Xw=         (Q0.*S./(m*U0)).*(CL_0-CDalp);

Mw=             (Q0.*S*c./(U0.*Iyy))*(Cmalp);

Zw_dot=         -(Q0.*S*c./(2*m.*U0.^2))*(CLalp_dot);

Mw_dot=         (Q0.*S*c^2./(2.*Iyy.*U0.^2))*(Cmalp_dot);

Zq=             -(Q0.*S*c./(2*m.*U0))*(CLq);

Mq=         (Q0.*S*c^2./(2.*Iyy.*U0))*(Cmq);

end
