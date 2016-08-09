import numpy as np
import matplotlib.pyplot as plt
import scipy.special as S
         
data1='/Users/LJ/output_gas_vphi.txt'
gas=[]   
f=open(data1,'r')
f.seek(0)
for i in range(4194304):
 List=f.readline().strip().split()
 gas.append(List)

data2='/Users/LJ/output_interpolation.txt'
star=[]
f=open(data2,'r')
f.seek(0)
for k in range(4194304):
 List=f.readline().strip().split()
 star.append(List)

data3='/Users/LJ/tangent_point.txt'
point=[]
f=open(data3,'r')
f.seek(0)
for k in range(900):
 List=f.readline().strip().split()
 point.append(List)

v_star=[]
v_gas=[]
for i in range(4194304):
 v_star.append(float(star[i][2]))
 v_gas.append(float(gas[i][2]))

v_star=np.asarray(v_star)
v_gas=np.asarray(v_gas)
diff=-v_gas-v_star

difference=diff.reshape((2048,2048))
difference=np.transpose(difference)

y=[]
for i in range(2048):
 y.append(gas[i][1])

x=y

x_plot1=[]
x_plot2=[]
y_plot1=[]
y_plot2=[]
for i in range(900):
 if 0.67<float(point[i][3])<=3.5:
  x_plot1.append(point[i][0])
  y_plot1.append(point[i][1])
 if 0.47<float(point[i][7])<=2.6:
  x_plot2.append(point[i][4])
  y_plot2.append(point[i][5])

plt.figure()
CS=plt.contourf(x,y,difference,50)
cbar=plt.colorbar(CS)
plt.plot(x_plot1,y_plot1,'.w')
plt.plot(x_plot2,y_plot2,'.w')
plt.title("phi Velocity Difference between Gas and Star")
plt.xlabel("x[kpc]")
plt.ylabel("y[kpc]")
plt.show()
