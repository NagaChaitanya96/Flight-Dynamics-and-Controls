function dLpdt=LP_fsfb(t,Lp)
global Alp_w_control a11 a12 a21 a22 ef
dLpdt=zeros(2,1);
u_cl=Lp(1);theta_cl=Lp(2);

tlp=1; Dtlp=20;
if ef==2
    if t>=tlp && t<=tlp+Dtlp
        u_cl=u_cl+0;
        theta_cl=theta_cl+0.087266; % giving disturbance 5degrees in radians 
    end
end
dLpdt(1)=a11*u_cl+a12*theta_cl;
dLpdt(2)= a21*u_cl+ a22*theta_cl;
end