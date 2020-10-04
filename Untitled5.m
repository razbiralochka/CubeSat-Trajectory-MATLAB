% By Khairullin Ilnar and Samadov Mahdi
clc;
clear all;
j=1;
for i=300:10:760
    ro(j,1)=(i);
    j=j+1;
end;
ro (1,2)=3.469/10^2;
ro (2,2)=2.854/10^2;
ro (3,2)=2.360/10^2;
ro (4,2)=1.960/10^2;
ro (5,2)=1.635/10^2;
ro (6,2)=1.369/10^2;
ro (7,2)=1.151/10^2;
ro (8,2)=9.704/10^3;
ro (9,2)=8.212/10^3;
ro (10,2)=6.970/10^3;
ro (11,2)=5.934/10^3;
ro (12,2)=5.066/10^3;
ro (13,2)=4.337/10^3;
ro (14,2)=3.722/10^3;
ro (15,2)=3.201/10^3;
ro (16,2)=2.760/10^3;
ro (17,2)=2.365/10^3;
ro (18,2)=2.065/10^3;
ro (19,2)=1.702/10^3;
ro (20,2)=1.557/10^3;
ro (21,2)=1.356/10^3;
ro (22,2)=1.183/10^3;
ro (23,2)=1.034/10^3;
ro (24,2)=0.053/10^4;
ro (25,2)=7.037/10^3;
ro (26,2)=6.070/10^4;
ro (27,2)=6.120/10^4;
ro (28,2)=5.308/10^4;
ro (29,2)=4.760/10^4;
ro (30,2)=4.204/10^4;
ro (31,2)=3.717/10^4;
ro (32,2)=3.377/10^4;
ro (33,2)=3.044/10^4;
ro (34,2)=2.747/10^4;
ro (35,2)=2.482/10^4;
ro (36,2)=2.246/10^4;
ro (37,2)=2.035/10^4;
ro (38,2)=1.846/10^4;
ro (39,2)=1.676/10^4;
ro (40,2)=1.524/10^4;
ro (41,2)=1.387/10^4;
ro (42,2)=1.583/10^4;
ro (43,2)=1.152/10^4;
ro (44,2)=1.051/10^4;
ro (45,2)=0.605/10^5;
ro (46,2)=8.783/10^5;
ro (47,2)=8.038/10^5;
S0=30*(pi/180); %star time on the Greenwich Meridian at Greenwich midnight, rad
Cx=2.2; %drag coefficient
gamma_min=10*(pi/180); %the minimum elevation angle, rad
ts0=0;
w=7.292115*10^(-5); %angular velocity of the Earth's rotation, 1/s
m=4; %the mass of SPACECRAFT, kg
Unit=3; %number of units, PCs
my=398602; %gravitational parameter of the Earth, km3/S2
Rz=6378.1; %the radius of Earth, km
Hp=700; %perigee height, km
Ha=700; %apogee height, km                  
Rp=Rz+Hp; %perigee radius, km
Ra=Rz+Ha; %apogee radius, km
i=45*(pi/180);% orbital inclination, rad
Omega=10*(pi/180); % longitude of the ascending node, rad
omega=0; % argument of perigee, rad
Ia=0; % true anomaly, rad
e=(Ra-Rp)/(Ra+Rp); %eccentricity
p=2*Ra*Rp/(Ra+Rp); %focal parameter
u=omega+Ia; %argument of latitude, rad
Smax=0.0001*sqrt(0.0001^2+0.0003^2);
sigma=Cx*Smax/m;
r=p/(1+e*cos(Ia)); %radius-vector
Vr=sqrt(my/p)*e*sin(Ia); %the radial component of velocity
Vn=sqrt(my/p)*(1+e*cos(Ia)); %transversal component of velocity
x0=r*(cos(u)*cos(Omega)-sin(u)*cos(i)*sin(Omega));
y0=r*(cos(u)*sin(Omega)+sin(u)*cos(i)*cos(Omega));
z0=r*sin(u)*sin(i);
vx0=Vr*(cos(u)*cos(Omega)-sin(u)*cos(i)*sin(Omega))-Vn*(cos(Omega)*sin(u)+sin(Omega)*cos(u)*cos(i));
vy0=Vr*(cos(u)*sin(Omega)+sin(u)*cos(i)*cos(Omega))-Vn*(sin(Omega)*sin(u)-cos(Omega)*cos(u)*cos(i));
vz0=Vr*sin(u)*sin(i)+Vn*cos(u)*sin(i);
rv0=[];
rv0(1)=x0;
rv0(2)=y0;
rv0(3)=z0;
rv0(4)=vx0;
rv0(5)=vy0;
rv0(6)=vz0;
%Runge-Kutt 4 orders
t0=0;
tk=100000;
h=50;
 [T,RV,RVV]=Rangikut4( h,t0,rv0,tk,sigma,ro );
%Calculation of the radius-vector and velocity modulus at each time point
rv=sqrt(RVV(:,1).^2+RVV(:,2).^2+RVV(:,3).^2);
Vv=sqrt(RVV(:,4).^2+RVV(:,5).^2+RVV(:,6).^2);
r=sqrt(RV(:,1).^2+RV(:,2).^2+RV(:,3).^2);
V=sqrt(RV(:,4).^2+RV(:,5).^2+RV(:,6).^2);
%Calculation of the true anomaly at each time point
n=size(RVV);
n=n(1,1);
rv_gr=[];
%Transfer to Greenwich UK
for k=1:n
    S=S0+w*(T(k,1)-ts0);
    Xg=RVV(k,1)*cos(S)+RVV(k,2)*sin(S);
    Yg=-RVV(k,1)*sin(S)+RVV(k,2)*cos(S);
    Zg=RVV(k,3);
    Vxg=RVV(k,4)*cos(S)+RVV(k,5)*sin(S)+w*RVV(k,2);
    Vyg=-RVV(k,4)*sin(S)+RVV(k,5)*cos(S)-w*RVV(k,1);
    Vzg=RVV(k,6);
    rv_gr=[rv_gr; Xg Yg Zg Vxg Vyg Vzg];
    fi(k)=atand(Zg/sqrt(Xg^2+Yg^2));%latitude, deg, [-i;i]
    la(k)=atan2d(Yg,Xg);%longitude, deg, [-180;180] 
end;
r_gr=sqrt(rv_gr(:,1).^2+rv_gr(:,2).^2+rv_gr(:,3).^2);
%Calculation of the visibility zone
fi_np=53.2001*(pi/180);%latitude of the ground point, rad
la_np=50.15*(pi/180); %longitude of the ground point, rad
X_np=Rz*cos(fi_np)*cos(la_np);
Y_np=Rz*cos(fi_np)*sin(la_np);
Z_np=Rz*sin(fi_np);
Vremya=0;
for k=1:n
    F(k,1)=(rv_gr(k,1)-X_np)*X_np+(rv_gr(k,2)-Y_np)*Y_np+(rv_gr(k,3)-Z_np)*Z_np-sqrt((rv_gr(k,1)-X_np)^2+(rv_gr(k,2)-Y_np)^2+(rv_gr(k,3)-Z_np)^2)*Rz*sin(gamma_min);
    if F>0
       Vremya=Vremya+h;
    end;
end;
% EarthImage
original=imread('Earth4.jpg');

N = 30;
phivec = linspace(0,pi,N);
thetavec = linspace(0,2*pi,2*N);
[ph, th] = meshgrid(thetavec,phivec);
R=6378.1*ones(size(th)); % should be your R(theta,phi) surface in general
x = R.*sin(th).*cos(ph);
y = R.*sin(th).*sin(ph);
z = R.*cos(th);
figure;
surf(x,y,z);
axis vis3d
%%%
h=findobj('Type','surface');
set(h,'CData',original,'FaceColor','texturemap')
hold on;
% % % 
%Orbit of the spacecraft
plot3(RVV(:,1), RVV(:,2), RVV(:,3),"red");
'Орбита космического аппрата';
xlabel('X, км'); 
ylabel('Y, км');
zlabel('Z, км');
grid on;
% % % 
%Plotting
%Dependence of the radius-vector from time to time
figure();
plot(T,r-rv, "red");
xlim([0 30000]);
xlabel('Время-t, секунды'); 
ylabel('Разница радиусов-векторов r-rv, км'); 
grid on;
%KA route
figure();
imagesc([-180 180], [-80 80], original);
hold on;
plot(la(1,1:601),0-fi(1,1:601), "red");
xlabel('Долгота-la, градус'); 
ylabel('Широта-fi, градус'); 
grid on;






