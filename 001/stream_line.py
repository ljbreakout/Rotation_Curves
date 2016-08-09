import numpy as np
import matplotlib.pyplot as plt

data1='/Users/LJ/output_gas_vphi.txt'
v_phi=[]
f=open(data1,'r')
f.seek(0)
for i in range(4194304):
 List=f.readline().strip().split()
 v_phi.append(List)

vphi=[]
for i in range(4194304):
 vphi.append(-float(v_phi[i][2]))

y_vphi=[]
for i in range(2048):
 y_vphi.append(float(v_phi[i][1]))

x_vphi=y_vphi
x_vphi=np.asarray(x_vphi)
y_vphi=np.asarray(y_vphi)
vphi=np.asarray(vphi)
vel_phi=vphi.reshape((2048,2048))
vel_phi=np.transpose(vel_phi)

data2='/Users/LJ/gas_kinematicsdat.txt'
gas=[]
f=open(data2,'r')
f.seek(0)
for i in range(4194304):
 List=f.readline().strip().split()
 gas.append(List)

y_gas=y_vphi
print y_gas
print type(y_gas)
x_gas=y_gas
vx_ori=[]
vy_ori=[]
for i in range(4194304):
 vx_ori.append(float(gas[i][3]))
 vy_ori.append(float(gas[i][4]))
vx_ori=np.asarray(vx_ori)
vy_ori=np.asarray(vy_ori)
vx_gas=vx_ori.reshape((2048,2048))
vx_gas=np.transpose(vx_gas)
vy_gas=vy_ori.reshape((2048,2048))
vy_gas=np.transpose(vy_gas)

x1=[]
y1=[]
x=-1.56006
y=-1.01807
vx=220.049
vy=-90.0326
t=1.5*10**12
for i in range(10000):
 x=x+float(vx)/(3.08*10**16)*t
 x1.append(x)
 y=y+float(vy)/(3.08*10**16)*t
 y1.append(y)
 index1=np.where(x_gas<(x))
 index2=np.where(x_gas==max(x_gas[index1]))
 index3=np.where(x_gas>=x)
 index4=np.where(x_gas==min(x_gas[index3]))
 index5=np.where(y_gas<y)
 index6=np.where(y_gas==max(y_gas[index5]))
 index7=np.where(y_gas>=y)
 index8=np.where(y_gas==min(y_gas[index7]))
 vx11=vx_gas[index6,index2]
 vx12=vx_gas[index8,index2]
 vx21=vx_gas[index6,index4]
 vx22=vx_gas[index8,index4]
 coef0=1/((x_gas[index4]-x_gas[index2])*(y_gas[index8]-y_gas[index6]))
 coef1=vx11*(x_gas[index4]-x)*(y_gas[index8]-y)
 coef2=vx21*(x-x_gas[index2])*(y_gas[index8]-y)
 coef3=vx12*(x_gas[index4]-x)*(y-y_gas[index6])
 coef4=vx22*(x-x_gas[index2])*(y-y_gas[index6])
 vx=coef0*(coef1+coef2+coef3+coef4)

 vy11=vy_gas[index6,index2]
 vy12=vy_gas[index8,index2]
 vy21=vy_gas[index6,index4]
 vy22=vy_gas[index8,index4]
 coef11=vy11*(x_gas[index4]-x)*(y_gas[index8]-y)
 coef22=vy21*(x-x_gas[index2])*(y_gas[index8]-y)
 coef33=vy12*(x_gas[index4]-x)*(y-y_gas[index6])
 coef44=vy22*(x-x_gas[index2])*(y-y_gas[index6])
 vy=coef0*(coef11+coef22+coef33+coef44)
 print i
plt.figure()
CS=plt.contourf(x_vphi,y_vphi,vel_phi,30)
cbar=plt.colorbar(CS)
plt.plot(x1,y1,'k')
plt.xlim(-5.0,5.0)
plt.ylim(-5.0,5.0)
plt.title("Stream Line")
plt.xlabel("x[kpc]")
plt.ylabel("y[kpc]")
plt.legend()
plt.show()
