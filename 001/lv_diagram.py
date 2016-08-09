import numpy as np
import matplotlib.pyplot as plt
import math

data1='/Users/LJ/gas_kinematicsdat.txt'
gas=[]  
f=open(data1,'r')
f.seek(0)    
for i in range(4194304):
 List=f.readline().strip().split()
 gas.append(List)
    
x_sun=8.3*np.cos(153./180.*np.pi)
y_sun=8.3*np.sin(153./180.*np.pi)
v_sun=210.
v_sun_x=-v_sun*np.sin(27./180.*np.pi)
v_sun_y=-v_sun*np.cos(27./180.*np.pi)
l=[]

for i in range(4194304):
 if np.sqrt((float(gas[i][0])**2+float(gas[i][1])**2))>=14.5:
  c=0
  l.append(c)
 if np.sqrt((float(gas[i][0])**2+float(gas[i][1])**2))<14.5:
  if float(gas[i][1])>(y_sun/x_sun)*float(gas[i][0]) and abs(float(gas[i][1])-(y_sun/x_sun)*float(gas[i][0]))>10**(-6):
   c=-1*math.acos((x_sun**2+y_sun**2-x_sun*float(gas[i][0])-y_sun*float(gas[i][1]))/(np.sqrt(x_sun**2+y_sun**2)*np.sqrt((float(gas[i][0])-x_sun)**2+(float(gas[i][1])-y_sun)**2)))*180/np.pi
   l.append(c)
  if float(gas[i][1])<(y_sun/x_sun)*float(gas[i][0]) and abs(float(gas[i][1])-(y_sun/x_sun)*float(gas[i][0]))>10**(-6):
   c=1*math.acos((x_sun**2+y_sun**2-x_sun*float(gas[i][0])-y_sun*float(gas[i][1]))/(np.sqrt(x_sun**2+y_sun**2)*np.sqrt((float(gas[i][0])-x_sun)**2+(float(gas[i][1])-y_sun)**2)))*180/np.pi
   l.append(c)
  if abs(float(gas[i][1])-(y_sun/x_sun)*float(gas[i][0]))<=10**(-6):
   c=0
   l.append(c)

v_los=[]
omega_b=33.
tmp_r=np.zeros(4194304)
e_los_x=np.zeros(4194304)
e_los_y=np.zeros(4194304)
velocity_x=np.zeros(4194304)
velocity_y=np.zeros(4194304)
v_los_tmp=np.zeros(4194304)
v_los_sun=np.zeros(4194304)

for k in range(4194304):
 if np.sqrt((float(gas[k][0])**2+float(gas[k][1])**2))<14.5:
  tmp_r[k]=np.sqrt((float(gas[k][0])-x_sun)**2+(float(gas[k][1])-y_sun)**2)
  e_los_x[k]=(float(gas[k][0])-x_sun)/tmp_r[k]
  e_los_y[k]=(float(gas[k][1])-y_sun)/tmp_r[k]
  v_los_tmp[k]=float(gas[k][3])*e_los_x[k]+float(gas[k][4])*e_los_y[k]
  v_los_sun[k]=v_sun_x*e_los_x[k]+v_sun_y*e_los_y[k]
  c=v_los_tmp[k]-v_los_sun[k]
  v_los.append(c)
 if np.sqrt((float(gas[k][0])**2+float(gas[k][1])**2))>=14.5:
  c=0
  v_los.append(c)

l=np.asarray(l)
v_los=np.asarray(v_los)
l1=[]
v_los1=[]
for j in range(900):
 index1=np.where((l>=-90.0+0.1*j) & (l<-90.0+0.1*(j+1)))
 index2=np.where(v_los==min(v_los[index1]))
 l1.append(l[index2])
 v_los1.append(v_los[index2])

l2=[]
v_los2=[]
v2=[]
for n in range(900):
 index3=np.where((l>=0.0+0.1*n) & (l<0.0+0.1*(n+1)))
 index4=np.where(v_los==max(v_los[index3]))
 l2.append(l[index4])
 v_los2.append(v_los[index4])

fig=plt.figure()
plt.plot(l1,v_los1,'.')
plt.plot(l2,v_los2,'.')
plt.xlabel("l(degree)")
plt.ylabel("v_los(km/s)")
plt.title("lv diagram")
plt.xlim(-90.0,90.0)
plt.show()
