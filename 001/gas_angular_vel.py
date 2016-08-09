import numpy as np
import matplotlib.pyplot as plt
import math

data1='/Users/LJ/output_gas_vphi.txt'
gas=[]  
f=open(data1,'r')
f.seek(0)    
for i in range(4194304):
 List=f.readline().strip().split()
 gas.append(List)

x=[]
y=[]
vphi=[]
for i in range(4194304):
 x.append(float(gas[i][0]))
 y.append(float(gas[i][1]))
 vphi.append(-float(gas[i][2]))

x=np.asarray(x)
y=np.asarray(y)
vphi=np.asarray(vphi)

R=np.sqrt(x*x+y*y)
omega1=vphi/R
omega=[]
for i in range(4194304):
 c=math.log10(omega1[i])
 omega.append(c)

omega=np.asarray(omega)

y_gas=[]
for i in range(2048):
 y_gas.append(y[i])

x_gas=y_gas
omega_gas=omega.reshape((2048,2048))
omega_gas=np.transpose(omega_gas)

#index=np.where((y==0.007324) & (R>0.1))
#R_plot=R[index]
#omega_plot=omega[index]

fig=plt.figure()
CS=plt.contourf(x_gas,y_gas,omega_gas,30)
cbar=plt.colorbar(CS)
#plt.plot(R_plot,omega_plot)
plt.xlabel("x[kpc]")
plt.ylabel("y[kpc]")
plt.xlim(-7.0,7.0)
plt.ylim(-7.0,7.0)
plt.title("Angular Velocity of Gas")
plt.show()
