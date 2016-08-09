import numpy as np
import scipy.special as S
import matplotlib.pyplot as plt 
         
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
 v1=np.sqrt(a1*R_star1)
 v_star.append(v1)

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

y=[]
for i in range(1500):
 y.append(star[i][1])

x=y

x=np.asarray(x)
y=np.asarray(y)
#v_star=np.asarray(v_star)
vt=np.asarray(vt)

vel_star=vt.reshape((1500,1500))
vel_star=np.transpose(vel_star)

plt.figure()
CS=plt.contourf(x,y,vel_star,30)
cbar=plt.colorbar(CS)
plt.title("phi Velocity of Star")
plt.xlabel("x[kpc]")
plt.ylabel("y[kpc]")
plt.show()
