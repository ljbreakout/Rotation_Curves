import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.special as S

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

#R_mod=[]
v_mod=[]
for i in range(4194304):
# x22=float(gas[i][0])**2
# y22=float(gas[i][1])**2
# r_mod=np.sqrt(x22+y22)
# R_mod.append(r_mod)
 vx22=float(gas[i][3])**2
 vy22=float(gas[i][4])**2
 v_mod1=np.sqrt(vx22+vy22)
 v_mod.append(v_mod1)
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

#print l

v_los=[]
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
#fig=plt.figure()
#plt.plot(l,v_los,'.')
#plt.xlim(-90.0,90.0)
#plt.show()

l=np.asarray(l)
v_los=np.asarray(v_los)
#R_mod=np.asarray(R_mod)
v_mod=np.asarray(v_mod)

R1=[]
R2=[]
v1=[]
v2=[]
v_mod1=[]

for j in range(900):
 index1=np.where((l>=-90.0+0.1*j) & (l<-90.0+0.1*(j+1)))
 index2=np.where(v_los==min(v_los[index1]))
 R=abs(8.3*np.sin(l[index2]/180.*np.pi))
 R1.append(R)
 v_mod1.append(v_mod[index2])
 v=v_los[index2]+v_sun*np.sin(l[index2]*np.pi/180)
 v1.append(-v)

fig=plt.figure()
plt.plot(R1,v1,'.r',label="$TP$")
plt.plot(R1,v_mod1,'.g',label="$model$",)
plt.xlabel("Radius(kpc)")
plt.ylabel("Velocity(km/s)")
#plt.xlim(0.0,1.7)
plt.title("V Difference between TP and Model(-90~0)")
plt.legend()
plt.show()

R2=[]
v2=[]
v_mod2=[]
for n in range(900):
 index3=np.where((l>=0.0+0.1*n) & (l<0.0+0.1*(n+1)))
 index4=np.where(v_los==max(v_los[index3]))
 R=abs(8.3*np.sin(l[index4]/180.*np.pi))
 R2.append(R)
 v=v_los[index4]+v_sun*np.sin(l[index4]*np.pi/180)
 v2.append(v)
 v_mod2.append(v_mod[index4])

#fig=plt.figure()
#plt.plot(R2,v2,'.r',label="$TP$")
#plt.plot(R2,v_mod2,'.g',label="$model$",)
#plt.xlabel("Radius(kpc)")
#plt.ylabel("Velocity(km/s)")
#plt.title("V Difference between TP and Model(0~90)")
#plt.xlim(0.0,1.7)
#plt.legend()
#plt.show()
