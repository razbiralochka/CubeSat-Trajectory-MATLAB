function [rv]= helpv(rv0,sigma,ro)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
mu=398602; 
Rz=6378.1;
C20=-0.001082626;
C40=2.371/1000000;
%C20=0;
%C40=0;
V=sqrt(rv0(4)^2+rv0(5)^2+rv0(6)^2);
r=sqrt(rv0(1)^2+rv0(2)^2+rv0(3)^2);
H=r-Rz;
i=1;
while ro(i,1)<H
    i=i+1;
end;
Ostatok=H-ro(i,1);
ro0=ro(i,2)+(ro(i+1,2)-ro(i,2))*Ostatok/10;
    rv(1)=rv0(4);
    rv(2)=rv0(5);
    rv(3)=rv0(6);
    rv(4)=-mu*rv0(1)/(r^3)+((3*C20*mu*(Rz^2)/((r^4)*2))*(5*(rv0(3)^2)/(r^2)-1)+(5*C40*mu*(Rz^4)/(8*(r^6)))*(3+(63*(rv0(3)^2)/(r^2)-42)*(rv0(3)^2)/(r^2)))*rv0(1)/r-sigma*ro0*V*rv0(4);
    rv(5)=-mu*rv0(2)/(r^3)+((3*C20*mu*(Rz^2)/((r^4)*2))*(5*(rv0(3)^2)/(r^2)-1)+(5*C40*mu*(Rz^4)/(8*(r^6)))*(3+(63*(rv0(3)^2)/(r^2)-42)*(rv0(3)^2)/(r^2)))*rv0(2)/r-sigma*ro0*V*rv0(5);
    rv(6)=-mu*rv0(3)/(r^3)+((3*C20*mu*(Rz^2)/((r^4)*2))*(5*(rv0(3)^2)/(r^2)-3)+(5*C40*mu*(Rz^4)/(8*(r^6)))*(15+(63*(rv0(3)^2)/(r^2)-70)*(rv0(3)^2)/(r^2)))*rv0(3)/r-sigma*ro0*V*rv0(6);
end
