
# 3D Heatmap in Python using matplotlib
  
# to make plot interactive 
  
# importing required libraries
from math import factorial
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from pylab import *

def radialN10State(r,n):
    return (2/(n*a_0))**(3/2)*(factorial(n-2)*factorial(n+1)/2/n)**(1/2)*np.exp(-r/(n*a_0))*(2*r/(n*a_0))*sum([(-2*r/(n*a_0))**i/(factorial(i)*factorial(n-2-i)*factorial(3+i))for i in range(n-1)])
def azimuthalN10State(theta):
    return 6**(1/2)/2*cos(theta)
def rThetaFromXYZ(x,y,z):
    return np.sqrt(x**2+y**2+z**2),np.arctan2(np.sqrt(x**2+y**2),z)
def totalN10State(x,y,z,n):
    return 1/np.sqrt(2*np.pi)*azimuthalN10State(np.arctan2(np.sqrt(x**2+y**2),z))*radialN10State(np.sqrt(x**2+y**2+z**2),n)
def total100State(x,y,z):
    return 1/(np.sqrt(np.pi)*a_0**(3/2))*np.exp(-np.sqrt(x**2+y**2+z**2)/a_0)
# creating a dummy dataset
nx,ny,nz=(10,10,10)
m_e = 9.109e-31
h_bar=1.005e-34
e=1.6e-19
epsilon_0=8.854e-12
a_0 = 4*np.pi*epsilon_0*h_bar**2/(m_e*e**2)
print(a_0)
cubeLeng = 2 # this is how many radii out to go (a_0)
x_min=-cubeLeng*a_0
x_max=cubeLeng*a_0
y_min=-cubeLeng*a_0
y_max=cubeLeng*a_0
z_min=-cubeLeng*a_0
z_max=cubeLeng*a_0
n = 29
nx = n
ny = n
nz = n
x = np.linspace(x_min,x_max,nx)
y = np.linspace(y_min,y_max,ny)
z = np.linspace(z_min,z_max,nz)
xx,yy,zz=np.meshgrid(x,y,z)
# creating figures
fig = plt.figure(figsize=(10, 10))
xx=xx.flatten()
yy=yy.flatten()
zz=zz.flatten()
alpha_value = 0.01
ax = fig.add_subplot(221, projection='3d')
colo = (total100State(xx,yy,zz))**2
colo = colo/np.max(colo)
color_map = cm.ScalarMappable(cmap=cm.hsv)
color_map.set_array(colo)
img = ax.scatter(xx, yy, zz, marker='s',s=10,c=cm.hsv(colo),alpha=alpha_value)
# plt.colorbar(color_map)
# ax.set_title("3D Heatmap")
ax.set_xticks(np.linspace(-cubeLeng*a_0,cubeLeng*a_0,5))
ax.set_xticklabels(map(str,range(-cubeLeng,cubeLeng+1,cubeLeng//2)))
ax.set_yticks(np.linspace(-cubeLeng*a_0,cubeLeng*a_0,5))
ax.set_yticklabels(map(str,range(-cubeLeng,cubeLeng+1,cubeLeng//2)))
ax.set_zticks(np.linspace(-cubeLeng*a_0,cubeLeng*a_0,5))
ax.set_zticklabels(map(str,range(-cubeLeng,cubeLeng+1,cubeLeng//2)))
ax.set_xlabel('X (a_0)')
ax.set_ylabel('Y (a_0)')
ax.set_zlabel('Z (a_0)')

ax = fig.add_subplot(222, projection='3d')
colo = (totalN10State(xx,yy,zz,2))**2
colo = colo/np.max(colo)
color_map = cm.ScalarMappable(cmap=cm.hsv)
color_map.set_array(colo)
img = ax.scatter(xx, yy, zz, marker='s',s=10,c=cm.hsv(colo),alpha=alpha_value)
# plt.colorbar(color_map)
# ax.set_title("3D Heatmap")
ax.set_xticks(np.linspace(-cubeLeng*a_0,cubeLeng*a_0,5))
ax.set_xticklabels(map(str,range(-cubeLeng,cubeLeng+1,cubeLeng//2)))
ax.set_yticks(np.linspace(-cubeLeng*a_0,cubeLeng*a_0,5))
ax.set_yticklabels(map(str,range(-cubeLeng,cubeLeng+1,cubeLeng//2)))
ax.set_zticks(np.linspace(-cubeLeng*a_0,cubeLeng*a_0,5))
ax.set_zticklabels(map(str,range(-cubeLeng,cubeLeng+1,cubeLeng//2)))
ax.set_xlabel('X (a_0)')
ax.set_ylabel('Y (a_0)')
ax.set_zlabel('Z (a_0)')

ax = fig.add_subplot(223, projection='3d')
colo = (totalN10State(xx,yy,zz,3))**2
colo = colo/np.max(colo)
color_map = cm.ScalarMappable(cmap=cm.hsv)
color_map.set_array(colo)
img = ax.scatter(xx, yy, zz, marker='s',s=10,c=cm.hsv(colo),alpha=alpha_value)
# plt.colorbar(color_map)
# ax.set_title("3D Heatmap")
ax.set_xticks(np.linspace(-cubeLeng*a_0,cubeLeng*a_0,5))
ax.set_xticklabels(map(str,range(-cubeLeng,cubeLeng+1,cubeLeng//2)))
ax.set_yticks(np.linspace(-cubeLeng*a_0,cubeLeng*a_0,5))
ax.set_yticklabels(map(str,range(-cubeLeng,cubeLeng+1,cubeLeng//2)))
ax.set_zticks(np.linspace(-cubeLeng*a_0,cubeLeng*a_0,5))
ax.set_zticklabels(map(str,range(-cubeLeng,cubeLeng+1,cubeLeng//2)))
ax.set_xlabel('X (a_0)')
ax.set_ylabel('Y (a_0)')
ax.set_zlabel('Z (a_0)')

ax = fig.add_subplot(224, projection='3d')
colo = (totalN10State(xx,yy,zz,4))**2
colo = colo/np.max(colo)
color_map = cm.ScalarMappable(cmap=cm.hsv)
color_map.set_array(colo)
img = ax.scatter(xx, yy, zz, marker='s',s=10,c=cm.hsv(colo),alpha=alpha_value)
# plt.colorbar(color_map)
# ax.set_title("3D Heatmap")
ax.set_xticks(np.linspace(-cubeLeng*a_0,cubeLeng*a_0,5))
ax.set_xticklabels(map(str,range(-cubeLeng,cubeLeng+1,cubeLeng//2)))
ax.set_yticks(np.linspace(-cubeLeng*a_0,cubeLeng*a_0,5))
ax.set_yticklabels(map(str,range(-cubeLeng,cubeLeng+1,cubeLeng//2)))
ax.set_zticks(np.linspace(-cubeLeng*a_0,cubeLeng*a_0,5))
ax.set_zticklabels(map(str,range(-cubeLeng,cubeLeng+1,cubeLeng//2)))
ax.set_xlabel('X (a_0)')
ax.set_ylabel('Y (a_0)')
ax.set_zlabel('Z (a_0)')

# displaying plot
plt.show()