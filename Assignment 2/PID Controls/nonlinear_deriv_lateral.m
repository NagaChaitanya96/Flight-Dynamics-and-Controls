function [Y_beta,Yp,Yr,N_beta,Np,Nr,L_beta,Lp,Lr,Y_deltaA,Y_deltaR,N_deltaA,N_deltaR,L_deltaA,L_deltaR,Yv,Lv,Nv,k,c1,c2,c3]=nonlinear_deriv_lateral()
global  Ixz Ixx Izz CY_beta CYp CYr CY_deltaR CY_deltaA Cn_beta Cl_beta Clp Clr Cl_deltaR Cl_deltaA   Cnp Cnr Cn_deltaR Cn_deltaA m Ixx Izz V b e S
rho=1.0581;
g=9.81;
AR=19.9;
Q=0.5.*rho.*V.^2;
k=1/(pi*AR*e);

%Lateral Derivatives
Y_beta=Q.*S.*CY_beta./m;
Yp=Q.*S.*b.*CYp./(2*m.*V);
Yr=Q.*S.*b.*CYr./(2*m.*V);

N_beta=Q.*S.*b.*Cn_beta./Izz;
Np=Q.*S.*b.^2.*Cnp./(2*Izz.*V);
Nr=Q.*S.*b.^2.*Cnr./(2*Izz.*V);

c1=Ixz/Ixx;
c2=Ixz/Izz;
c3=c1*c2;

L_beta=Q.*S.*b.*Cl_beta./Ixx;
Lp=Q.*S.*b.^2.*Clp./(2*Ixx.*V);
Lr=Q.*S.*b.^2.*Clr./(2*Ixx.*V);


Y_deltaA=Q.*S.*CY_deltaA./m;
Y_deltaR=Q.*S.*CY_deltaR./m;

N_deltaA=Q.*S.*b.*Cn_deltaA./Izz;
N_deltaR=Q.*S.*b.*Cn_deltaR./Izz;

L_deltaA=Q.*S.*b.*Cl_deltaA./Ixx;
L_deltaR=Q.*S.*b.*Cl_deltaR./Ixx;

Yv=Y_beta./V;
Lv=L_beta./V;
Nv=N_beta./V;
c1=Ixz/Ixx;
c2=Ixz/Izz;
c3=c1*c2;
end
