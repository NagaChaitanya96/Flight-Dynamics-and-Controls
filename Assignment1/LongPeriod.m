function [dLodt]=LongPeriod(t,Lo)
global Xu  Zu  U0 theta0 Zq ef j
g=9.81;
dLodt=zeros(2,1);
u=Lo(1); theta=Lo(2);
tlp=1; Dtlp=20;
if ef==2
    if t>=tlp && t<=tlp+Dtlp
        u=u+0;
        theta=theta+0.1;
    end
end
dLodt(1)=Xu(j)*u+(-g*cos(theta0))*theta;
dLodt(2)= (-Zu(j)./(U0(j)+Zq(j))).*u+ (g.*sin(theta0)./(U0(j)+Zq(j)))*theta;
end