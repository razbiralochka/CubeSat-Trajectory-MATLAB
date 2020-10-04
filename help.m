function [rv]= help(rv0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
r=sqrt(rv0(1)^2+rv0(2)^2+rv0(3)^2);
mu=398602; 
    rv(1)=rv0(4);
    rv(2)=rv0(5);
    rv(3)=rv0(6);
    rv(4)=-mu*rv0(1)/(r^3);
    rv(5)=-mu*rv0(2)/(r^3);
    rv(6)=-mu*rv0(3)/(r^3);
end
