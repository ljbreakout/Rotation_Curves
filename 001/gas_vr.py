import numpy as np
import matplotlib.pyplot as plt

data='/Users/LJ/gas_kinematicsdat.txt'
gas=[]   
f=open(data,'r')
f.seek(0)
for i in range(4194304):
 List=f.readline().strip().split()
 gas.append(List)

v_r=[]
for i in range(4194304): 
 R=np.sqrt(float(gas[i][0])**2+float(gas[i][1])**2)
 e_x=(float(gas[i][0]))/R
 e_y=(float(gas[i][1]))/R
 v_r1=float(gas[i][3])*(e_x)+float(gas[i][4])*(e_y)
 v_r.append(v_r1)

print min(v_r)
print max(v_r)
y=[] 
for i in range(2048):
 y.append(gas[i][1])

x=y
x=np.asarray(x)
y=np.asarray(y)
v_r=np.asarray(v_r)
vel_r=v_r.reshape((2048,2048))
vel_r=np.transpose(vel_r)

data3='/Users/LJ/tangent_point.txt'
point=[]
f=open(data3,'r')
f.seek(0)
for k in range(900):
 List=f.readline().strip().split()
 point.append(List)

x_plot1=[]
x_plot2=[]
y_plot1=[]
y_plot2=[]
for i in range(900):
 if np.sqrt(float(point[i][0])**2+float(point[i][1])**2)<=2.6:
  x_plot1.append(point[i][0])
  y_plot1.append(point[i][1])
 if np.sqrt(float(point[i][2])**2+float(point[i][3])**2)<=3.5:
  x_plot2.append(point[i][2])
  y_plot2.append(point[i][3])

plt.figure()
CS=plt.contourf(x,y,vel_r,50)
cbar=plt.colorbar(CS)
plt.plot(x_plot1,y_plot1,'.w')
plt.plot(x_plot2,y_plot2,'.w')
plt.xlim(-4.0,4.0)
plt.ylim(-4.0,4.0)
plt.title("R Velocity of Gas")
plt.xlabel("x[kpc]")
plt.ylabel("y[kpc]")
plt.show()
