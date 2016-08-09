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

R1=[]
R2=[]
v1=[]
v2=[]

for j in range(900):
 index1=np.where((l>=-90.0+0.1*j) & (l<-90.0+0.1*(j+1)))
 index2=np.where(v_los==min(v_los[index1]))
 R=abs(8.3*np.sin(l[index2]/180.*np.pi))
 R1.append(R)
 v=v_los[index2]+v_sun*np.sin(l[index2]*np.pi/180)
 v1.append(-v)

R1=np.asarray(R1)
v1=np.asarray(v1)

#for n in range(900):
# index3=np.where((l>=0.0+0.1*n) & (l<0.0+0.1*(n+1)))
# index4=np.where(v_los==max(v_los[index3]))
# R=abs(8.3*np.sin(l[index4]/180.*np.pi))
# R2.append(R)
# v=v_los[index4]+v_sun*np.sin(l[index4]*np.pi/180)
# v2.append(v)
#R2=np.asarray(R2)
#v2=np.asarray(v2)
data2='/Users/LJ/output_interpolation.txt'
star=[]
f=open(data2,'r')
f.seek(0)
for i in range(4194304):
 List=f.readline().strip().split()
 star.append(List)

R_star=[]
v_star=[]
for i in range(4194304):
 if ((float(star[i][1])==0.007324) & (float(star[i][0])>0)):
  R_star.append(float(star[i][0]))
  v_star.append(float(star[i][2]))
print len(R_star)
print len(v_star)
R_star=np.asarray(R_star)
v_star=np.asarray(v_star)
R_sele=[]
v_sele=[]
for i in range(900):
 diff=abs(R_star-R1[i])
 index1=np.where(diff==min(diff))
 R_sele.append(float(R_star[index1]))
 v_sele.append(float(v_star[index1]))

v_sele=np.asarray(v_sele)
modify=[]
R_plot=[]
for i in range(900):
 a=v_sele[i]/v1[i]
 modify.append(float(a))
 R_plot.append(float(R1[i]))

modify=np.asarray(modify)
R_plot=np.asarray(R_plot)
data=np.array([R_plot,modify])
data=data.T

myfile=open('/Users/LJ/output_modify1.txt','w+')
np.savetxt(myfile,data,fmt=['%f','%f'])
myfile.close()

fig=plt.figure()
plt.plot(R_plot,modify)
plt.xlabel("Radius(kpc)")
plt.ylabel("Difference")
plt.title("Modify Function(-90~0)")
plt.show()
