import numpy as np
import math as m
import matplotlib.pyplot as plt



DawnData=np.genfromtxt('DawnFull_tb.txt')*1000
EarthData=np.genfromtxt('Earth_tb.txt')*1000
MarsData=np.genfromtxt('Mars_tb.txt')*1000
VestaData=np.genfromtxt('Vesta_tb.txt')*1000
CeresData=np.genfromtxt('Ceres_tb.txt')*1000

##Dawn JPL coodrds
X=DawnData[:,0]
Y=DawnData[:,1]
Z=DawnData[:,2]

CalcX=[]
CalcY=[]
Thrust=np.zeros(2714)
LambdMas=np.zeros(2714)
days=np.arange(2714)+1



##Earth JPL coords
EarthX=EarthData[:,0]
EarthY=EarthData[:,1]
EarthZ=EarthData[:,2]

##Mars JPL Coords
MarsX=MarsData[:,0]
MarsY=MarsData[:,1]
MarsZ=MarsData[:,2]


##Vesta JPL Coords
VestaX=VestaData[:,0]
VestaY=VestaData[:,1]
VestaZ=VestaData[:,2]


##Ceres JPL Coords
CeresX=CeresData[:,0]
CeresY=CeresData[:,1]
CeresZ=CeresData[:,2]

SunParam=132712440018*pow(10,9);

radius = m.sqrt(X[0]*X[0]+Y[0]*Y[0])
radius0 = radius
print(radius0)
angle = m.atan(Y[0]/X[0])

V0=np.array([-3.267334027266680E+00,3.345804245220296E+01])*1000
mV0=m.sqrt(V0[0]*V0[0]+V0[1]*V0[1])
Vmars=np.array([2.239956000086868E+01, 1.417958296152212E+01])*1000

rV=(X[0]*V0[0]+Y[0]*V0[1])/radius

phiV=m.sqrt(mV0*mV0-rV*rV)

print(phiV, rV)


m0=747.1+425+45.6



c0=26000


fuel = 425
k = np.zeros(4)
q = np.zeros(4)
day = 1

dday = 1
mas=747.1+fuel+45.6






while day < 510:
    CalcX.append(radius*m.cos(angle))
    CalcY.append(radius*m.sin(angle))
    F0=3*50*0.000
    lambd=m.pi/2
    if day > 50 and day <300:
        F0=3*90*0.001 * pow((radius0/radius),2)     
        lambd=m.pi/2+m.pi/2.6
    if day > 299 and day <500:
        F0=3*90*0.001 * pow((radius0/radius),2)     
        lambd=m.pi/2
    
    
    
    LambdMas[day-1]=lambd
    Thrust[day-1]=F0
    a = F0/mas
    consumption=F0/c0
    
    angle = angle + (phiV/radius)*dday*24*3600
    radius = radius + rV * dday*24*3600
    rV=rV+(pow(phiV,2)/radius-SunParam/pow(radius,2)+m.cos(lambd)*a)*dday*24*3600
    phiV=phiV+(-(rV*phiV)/radius+m.sin(lambd)*a)*dday*24*3600
    mas=mas-consumption*dday*24*3600
    day = day + dday 
print(mas-747.1-45.6)   

radius = m.sqrt(MarsX[-1]*MarsX[-1]+MarsY[-1]*MarsY[-1])
angle = m.atan(MarsY[-1]/MarsX[-1])

rV=(MarsX[-1]*Vmars[0]+MarsY[-1]*Vmars[1])/radius
mV0=m.sqrt(Vmars[0]*Vmars[0]+Vmars[1]*Vmars[1])

phiV=m.sqrt(mV0*mV0-rV*rV)
print(phiV, rV)



##fly to Vesta

while day > 509 and day < 1357:
    CalcX.append(radius*m.cos(angle))
    CalcY.append(radius*m.sin(angle))
    F0=3*50*0.000
    lambd=m.pi/2
    if day > 723 and day <1250:
        F0=3*90*0.001 * pow((radius0/radius),2)      
        lambd=m.pi/2+m.pi/200
        
    
    
    
    LambdMas[day-1]=lambd
    Thrust[day-1]=F0
    a = F0/mas
    consumption=F0/c0
    
    angle = angle + (phiV/radius)*dday*24*3600
    radius = radius + rV * dday*24*3600
    rV=rV+(pow(phiV,2)/radius-SunParam/pow(radius,2)+m.cos(lambd)*a)*dday*24*3600
    phiV=phiV+(-(rV*phiV)/radius+m.sin(lambd)*a)*dday*24*3600
    mas=mas-consumption*dday*24*3600
    day = day + dday 
print(mas-747.1-45.6) 


##around the Vesta
while day > 1356 and day < 1804:
     lambd=m.pi/2
     LambdMas[day-1]=lambd
     CalcX.append(X[day])
     CalcY.append(Y[day])
     day = day + dday
print(mas-747.1-45.6)



##fly to Ceres
Vvesta=np.array([-1.430410330851633E+01, 1.056459146838846E+01])*1000
radius = m.sqrt(X[1804]*X[1804]+Y[1804]*Y[1804])
angle = m.atan(Y[1804]/X[1804])

rV=(X[1804]*Vvesta[0]+Y[1804]*Vvesta[1])/radius
mV0=m.sqrt(Vvesta[0]*Vvesta[0]+Vvesta[1]*Vvesta[1])

phiV=m.sqrt(mV0*mV0-rV*rV)
print(phiV, rV)



while day > 1803 and day < 2715:
    CalcX.append(radius*m.cos(angle))
    CalcY.append(radius*m.sin(angle))
    F0=3*50*0.000
    lambd=m.pi/2
    if day > 1890 and day <2300:
        F0=3*90*0.001 * pow((radius0/radius),2)      
        lambd=m.pi/2
   
    LambdMas[day-1]=lambd
    Thrust[day-1]=F0
    a = F0/mas
    consumption=F0/c0
    
    angle = angle + (phiV/radius)*dday*24*3600
    radius = radius + rV * dday*24*3600
    rV=rV+(pow(phiV,2)/radius-SunParam/pow(radius,2)+m.cos(lambd)*a)*dday*24*3600
    phiV=phiV+(-(rV*phiV)/radius+m.sin(lambd)*a)*dday*24*3600
    mas=mas-consumption*dday*24*3600
    day = day + dday 
print(mas-747.1-45.6) 
















#plt.plot(X,Y,color='black',label='Траектория КА Dawn JPL');

plt.plot(EarthX,EarthY,color='green',label='Орбита Земли')
plt.scatter(X[1],Y[1], color='green', s=40, marker='o')

plt.plot(CeresX,CeresY,color='blue',label='Орбита Цереры')
plt.scatter(CeresX[2715],CeresY[2715], color='blue', s=40, marker='o')

plt.plot(VestaX,VestaY,color='cyan',label='Орбита Весты')
plt.scatter(VestaX[-1],VestaY[-1], color='cyan', s=40, marker='o')
plt.scatter(X[1804],Y[1804], color='cyan', s=40, marker='x')
plt.plot(MarsX,MarsY,color='red',label='Орбита Марса')
plt.scatter(MarsX[-1],MarsY[-1], color='red', s=40, marker='o')

plt.plot(CalcX,CalcY,color='orange',label='Траектория КА Dawn Рассчёт');

plt.axis('equal')

plt.legend(loc=2)
##plt.grid()
plt.show()


fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(days,Thrust,'red')
ax2.plot(days,LambdMas*180/m.pi,'blue')
 
ax1.set_xlabel("X data")
ax1.set_ylabel("Тяга, мН",color='Red')
 
ax2.set_ylabel("угол поворота, гр",color='b')
plt.xlim([1, 2714])
plt.show()