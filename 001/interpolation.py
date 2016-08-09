import numpy as np
import matplotlib.pyplot as plt
import scipy.special as S

data='/Users/LJ/output_gas_vphi.txt'
gas=[]
f=open(data,'r')
f.seek(0)
for i in range(4194304):
 List=f.readline().strip().split()
 gas.append(List)

y_gas=[]
for i in range(2048):
 y_gas.append(float(gas[i][1]))

x_gas=y_gas
x_gas=np.asarray(x_gas)
y_gas=np.asarray(y_gas)

data2='/Users/LJ/output_star_vt.txt'
star=[]
f=open(data2,'r')
f.seek(0)
for k in range(2250000):
 List=f.readline().strip().split()
 star.append(List)

vt=[]
for h in range(2250000):
 vt.append(float(star[h][2]))

y_star=[]
for i in range(1500):
 y_star.append(float(star[i][1]))

x_star=y_star
x_star=np.asarray(x_star)
y_star=np.asarray(y_star)
vt=np.asarray(vt)
vel=vt.reshape((1500,1500))
vel_star=np.transpose(vel)

#calculate velocity phi of star

v_star_new=[]
o=0
for i in range(2048):
 index1=np.where(x_star<(x_gas[i]))
 index2=np.where(x_star==max(x_star[index1]))
 index3=np.where(x_star>=(x_gas[i]))
 index4=np.where(x_star==min(x_star[index3]))
 for j in range(2048):
  index5=np.where(y_star<(y_gas[j]))
  index6=np.where(y_star==max(y_star[index5]))
  index7=np.where(y_star>=(y_gas[j]))
  index8=np.where(y_star==min(y_star[index7]))
#pick up four points
  v11=vel_star[index6,index2]
  v12=vel_star[index8,index2]
  v21=vel_star[index6,index4]
  v22=vel_star[index8,index4]
#formula of bilinear interpolation
#pick up the vel value of four points
  coef0=1/((x_star[index4]-x_star[index2])*(y_star[index8]-y_star[index6]))
  coef1=v11*(x_star[index4]-x_gas[i])*(y_star[index8]-y_gas[j])
  coef2=v21*(x_gas[i]-x_star[index2])*(y_star[index8]-y_gas[j])
  coef3=v12*(x_star[index4]-x_gas[i])*(y_gas[j]-y_star[index6])
  coef4=v22*(x_gas[i]-x_star[index2])*(y_gas[j]-y_star[index6])
  v_new=coef0*(coef1+coef2+coef3+coef4)
  v_star_new.append(v_new)    
  o=o+1
  print o

print len(v_star_new)
x=[]
y=[]
for i in range(4194304):
 x.append(float(gas[i][0]))
 y.append(float(gas[i][1]))

x=np.asarray(x)
y=np.asarray(y)
v_star_new=np.asarray(v_star_new)
data=np.array([x,y,v_star_new])
data=data.T

myfile=open('/Users/LJ/output_interpolation.txt','w+')
np.savetxt(myfile,data,fmt=['%f','%f','%f'])
myfile.close()
