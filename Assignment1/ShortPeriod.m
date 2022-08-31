function dshdt= ShortPeriod(t,sh)
global Zw U0 Mw Zw_dot Mw_dot Zq Mq ef j
dshdt=zeros(2,1);
w= sh(1); q =sh(2);
tsp =1; Dtsp= 0.5;
if ef==1
    if t>=tsp && t<=tsp+Dtsp
        w=w+5;
    end
end
dshdt(1)= (Zw(j)./(1-Zw_dot(j)))*w+((U0(j)+Zq(j))./(1-Zw_dot(j)))*q;
dshdt(2)= (Mw(j)+(Mw_dot(j)*Zw(j))./(1-Zw_dot(j)))*w+ (Mq(j)+((U0(j)+Zq(j))*Mw_dot(j))./(1-Zw_dot(j)))*q;
end