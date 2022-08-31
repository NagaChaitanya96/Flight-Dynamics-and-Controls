function [dLo_ncdt]=Lp_without_c(t,Lo_nc)
global nc11 nc12 nc21 nc22 ef
dLo_ncdt=zeros(2,1);
u_ncl=Lo_nc(1);theta_ncl=Lo_nc(2);

tlp=1; Dtlp=20;
if ef==2
    if t>=tlp && t<=tlp+Dtlp
        u_ncl=u_ncl+0;
        theta_ncl=theta_ncl+0.087266; % giving disturbance 5degrees in radians 
    end
end
dLo_ncdt(1)=nc11*u_ncl+nc12*theta_ncl;
dLo_ncdt(2)= nc21*u_ncl+ nc22*theta_ncl;
end