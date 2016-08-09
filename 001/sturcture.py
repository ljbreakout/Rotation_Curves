import numpy as np
import matplotlib.pyplot as plt
import math

data='/Users/LJ/gas_kinematicsdat.txt'
gas=[]
f=open(data,'r')
f.seek(0)
for i in range(4194304):
 List=f.readline().strip().split()
 gas.append(List)

data1='/Users/LJ/tangent_point.txt'
TP=[]
f=open(data1,'r')
f.seek(0)
for i in range(900):
 List=f.readline().strip().split()
 TP.append(List)

x_sun=8.3*np.cos(153./180.*np.pi)
y_sun=8.3*np.sin(153./180.*np.pi)

y=[]
p=[]
for i in range(2048):
 y.append(gas[i][1])

x=y

for j in range(4194304):
 p1=math.log10(float(gas[j][2]))
 p.append(p1)

x=np.asarray(x)
y=np.asarray(y)
p=np.asarray(p)

den=p.reshape((2048,2048))
den=np.transpose(den)

x1=[]
y1=[]
R1=[]
x2=[]
y2=[]
R2=[]
for j in range(900):
 x1.append(float(TP[j][0]))
 y1.append(float(TP[j][1]))
 R1.append(float(TP[j][3]))
 x2.append(float(TP[j][4]))
 y2.append(float(TP[j][5]))
 R2.append(float(TP[j][7]))

x1=np.asarray(x1)
x2=np.asarray(x2)
R1=np.asarray(R1)
y1=np.asarray(y1)
y2=np.asarray(y2)
R2=np.asarray(R2)
index1=np.where((R1>4.68) & (R1<5.5))
x_stru1=x1[index1]
y_stru1=y1[index1]
index2=np.where((R2>4.68) & (R2<5.5))
x_stru2=x2[index2]
y_stru2=y2[index2]

x11=[x_sun,-7.243650]
y11=[y_sun,4.064940]
x22=[x_sun,-5.764160]
y22=[y_sun,4.826660]
x33=[x_sun,-0.065918]
y33=[y_sun,4.152830]
x44=[x_sun,0.710449]
y44=[y_sun,-0.373535]
x55=[x_sun,-4.753420]
y55=[y_sun,-0.314941]
x66=[x_sun,-7.185060]
y66=[y_sun,-0.270996]
x77=[x_sun,-7.624510]
y77=[y_sun,3.317870]
a1=np.polyfit(x11,y11,1)
z1=np.poly1d(a1)
plt.figure()
CS=plt.contourf(x,y,den,30)
cbar=plt.colorbar(CS)
plt.plot(x1,y1,'.w')
plt.plot(x2,y2,'.w')
plt.plot(x_stru1,y_stru1,'.k')
plt.plot(x_stru2,y_stru2,'.k')
plt.plot(x11,y11,'b')
plt.plot(x22,y22,'b')
plt.plot(x33,y33,'b')
plt.plot(x44,y44,'b')
plt.plot(x55,y55,'b')
plt.plot(x66,y66,'b')
plt.plot(x77,y77,'b')
plt.plot(x_sun,y_sun,'.k')
plt.xlim(-8.3,5)
plt.ylim(-5,8.3)
plt.title("Feature")
plt.xlabel("x[kpc]")
plt.ylabel("y[kpc]")
plt.legend()
plt.show()
