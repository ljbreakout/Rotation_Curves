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

#print l

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
#fig=plt.figure()
#plt.plot(l,v_los,'.')
#plt.xlim(-90.0,90.0)
#plt.show()

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
 #v=v_sun*np.sin(l[index2]*np.pi/180)-v_los[index2]
 v=v_los[index2]+v_sun*np.sin(l[index2]*np.pi/180)
 v1.append(-v)

#fig=plt.figure()
#plt.plot(R1,v1,'.')
#plt.show()

R2=[]
v2=[]
for n in range(900):
 index3=np.where((l>=0.0+0.1*n) & (l<0.0+0.1*(n+1)))
 index4=np.where(v_los==max(v_los[index3]))
 R=abs(8.3*np.sin(l[index4]/180.*np.pi))
 R2.append(R)
 v=v_los[index4]+v_sun*np.sin(l[index4]*np.pi/180)
 v2.append(v)

#fig=plt.figure()
#plt.plot(R2,v2,'.r')
#plt.show()

data2='/Users/LJ/potential_M80fit_15kpc.txt'
star=[]
f=open(data2,'r')
f.seek(0)
for k in range(2250000):
 List=f.readline().strip().split()
 star.append(List)

x_star=[]
y_star=[]
ax2=[]
ay2=[]
R_star=[]
v_star=[]

for h in range(2250000):
 x_star.append(star[h][0])
 y_star.append(star[h][1])
 ax2.append(star[h][6])
 ay2.append(star[h][7])

for h in range(2250000):
 R_star1=np.sqrt(float(x_star[h])*float(x_star[h])+float(y_star[h])*float(y_star[h]))
 R_star.append(R_star1)
 a1=np.sqrt(float(ax2[h])*float(ax2[h])+float(ay2[h])*float(ay2[h]))
 v11=np.sqrt(a1*R_star1)
 v_star.append(v11)


vt=[]

Mbh=3.5*10**7
eps=0.01
gconst=4.302*10**(-6)

par1=2.7856e-11
par2=1.34497e-12
coeff1=-3.65269e8
coeff2=-2.20803e8
const1=S.gamma(0.4)*(1.0-S.gammainc(0.4,(par1)*0.0))*coeff1
const2=S.gamma(0.4)*(1.0-S.gammainc(0.4,(par2)*0.0))*coeff2
mass2light=2.

for g in range(2250000):
 vbh2=gconst*Mbh*float(R_star[g])**2/((eps**2+float(R_star[g])**2)**1.5)
 tmp1_x=mass2light*0.5*(S.gamma(0.4)*(1.0-S.gammainc(0.4,(par1)*((float(R_star[g])*1000.)**5)))*(coeff1)-const1)
 tmp2_x=mass2light*0.5*(S.gamma(0.4)*(1.0-S.gammainc(0.4,(par2)*((float(R_star[g])*1000.)**5)))*(coeff2)-const2)
 vnb2=4.302*(tmp1_x+tmp2_x)/(float(R_star[g])*1.0e6)

 barmass=1.0*10**10
 a=9.0
 vlongbar2=gconst*barmass*float(R_star[g])/((float(R_star[g])+a)**2)
 vt1=np.sqrt((float(v_star[g]))**2+vbh2+vnb2+vlongbar2)
 vt.append(vt1)

R_plot=[]
vt_plot=[]
for q in range(2250000):
 if float(x_star[q])==0.1000595E-01:
  R_plot.append(float(R_star[q]))
  vt_plot.append(float(vt[q]))


fig=plt.figure()
plt.plot(R_plot,vt_plot,label="$true-RC$")
plt.plot(R1,v1,'.r',label="$gas--90-RC$")
plt.plot(R2,v2,'.g',label="$gas-90-RC$")
plt.xlabel("Radius(kpc)")
plt.ylabel("Velocity(km/s)")
plt.title("Rotation Curves")
plt.xlim(0.0,9.0)
plt.legend()
plt.show()
