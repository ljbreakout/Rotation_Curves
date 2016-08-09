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

v_phi=[]
p=[]
for i in range(4194304):
 R=np.sqrt(float(gas[i][0])**2+float(gas[i][1])**2)
 e_x=(float(gas[i][0]))/R
 e_y=(float(gas[i][1]))/R
 v_phi1=float(gas[i][3])*(e_y)+float(gas[i][4])*(-e_x)
 v_phi.append(v_phi1)
 p1=math.log10(float(gas[i][2]))
 p.append(p1)

print min(v_phi)
 
y=[]
for i in range(2048):
 y.append(gas[i][1])

x=y
x=np.asarray(x)
y=np.asarray(y)
v_phi=np.asarray(v_phi)

vel_phi=v_phi.reshape((2048,2048))
vel_phi=np.transpose(vel_phi)

p=np.asarray(p)
density=p.reshape((2048,2048))
density=np.transpose(density)

data2='/Users/LJ/tangent_point.txt'
point=[]
f=open(data2,'r')
f.seek(0)
for i in range(900):
 List=f.readline().strip().split()
 point.append(List)

x1=[]
y1=[]
x2=[]
y2=[]
for i in range(900):
 x1.append(float(point[i][0]))
 y1.append(float(point[i][1]))
 x2.append(float(point[i][4]))
 y2.append(float(point[i][5]))

plt.figure()
CS=plt.contourf(x,y,vel_phi,30)
cbar=plt.colorbar(CS)
plt.contour(x,y,density)
plt.plot(x1,y1,'.w')
plt.plot(x2,y2,'.w')
plt.title("phi Velocity of Gas")
plt.xlabel("x[kpc]")
plt.ylabel("y[kpc]")
plt.xlim(-8.3,5)
plt.ylim(-5,8.3)
plt.legend()
plt.show()
