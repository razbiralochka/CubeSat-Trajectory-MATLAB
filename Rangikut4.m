function [T,RV,RVv]=Rangikut4( h,t0,rv0,tk,sigma,ro )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    RVv=[];
    T(1,1)=t0;
    rv=rv0;
    RV=rv0;
    RVv=rv0;
    i=2;
    t=t0;
    while t<tk
        k1=helpv(rv0,sigma,ro);
        k2=helpv(rv0+h*k1/2,sigma,ro);
        k3=helpv(rv0+h*k2/2,sigma,ro);
        k4=helpv(rv0+h*k3,sigma,ro);
        k=(k1+2*k2+2*k3+k4)/6;
        rv0=rv0+h*k;
        t=t+h;
        RVv=[RVv; rv0];
        T(i,1)=t;
        i=i+1;
    end; 
    i=2;
    t=t0;
    while t<tk
        k1=help(rv);
        k2=help(rv+h*k1/2);
        k3=help(rv+h*k2/2);
        k4=help(rv+h*k3);
        k=(k1+2*k2+2*k3+k4)/6;
        rv=rv+h*k;
        t=t+h;
        RV=[RV; rv];
        i=i+1;
    end; 
 end
