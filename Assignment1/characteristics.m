function [dydt] = characteristics(t,y)
global  Xu  Zu  U0 theta0 Mu Zw Xw Mw Zw_dot Mw_dot Zq Mq  ef j
dydt=zeros(4,1);
g=9.81;
u= y(1); w=y(2) ; q=y(3); theta=y(4);
if ef==1
    tsp=1; Dtsp=0.5;
    if t>=tsp && t<= tsp+Dtsp
        w=w+5;
    end
elseif ef==2
    tlp=1; Dtlp=20;
    if t>=tlp && t<= tlp+Dtlp
        u=u+0;
        theta=theta+0.1;
    end
end
dydt(1)=Xu(j)*u+ Xw(j)*w + 0*q+ (-g*cos(theta0))*theta;
dydt(2)=(Zu(j)./(1-Zw_dot(j)))*u+ (Zw(j)./(1-Zw_dot(j)))*w+ ((U0(j)+Zq(j))./(1-Zw_dot(j)))*q -(g*sin(theta0)./(1-Zw_dot(j)))*theta;
dydt(3)=(Mu(j)+(Mw_dot(j)*Zu(j))./(1-Zw_dot(j)))*u+( Mw(j)+(Mw_dot(j)*Zw(j))./(1-Zw_dot(j)))*w+ (Mq(j)+((U0(j)+Zq(j))*Mw_dot(j))./(1-Zw_dot(j)))*q-(g*sin(theta0)*Mw_dot(j)./(1-Zw_dot(j)))*theta;
dydt(4)=q;
end