
# 3D Heatmap in Python using matplotlib
  
# to make plot interactive 
  
# importing required libraries
from math import factorial
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from math import comb
from matplotlib.widgets import Slider, Button
from sympy import E1

m_e = 9.109e-31
h_bar=1.055e-34
e=1.6e-19
epsilon_0=8.854e-12
a_0 = 4*np.pi*epsilon_0*h_bar**2/(m_e*e**2)
e_1=-h_bar**2/(2*m_e*a_0**2)
print(e_1)
def radialN10State(r,n):
    return (2/(n*a_0))**(3/2)*((n+1)*(n-1)/2)**(1/2)*np.exp(-r/(n*a_0))*(2*r/(n*a_0))*np.sum([comb(n-2,i)*(-2*r/(n*a_0))**i/factorial(3+i)for i in range(n-1)],axis=0)
def azimuthalN10State(theta):
    return 6**(1/2)/2*cos(theta)
def radial100State(r):
    return 2/(a_0**(3/2))*np.exp(-r/a_0)
def azimuthal100State(theta):
    return 1/np.sqrt(2)
def rThetaFromXYZ(x,y,z):
    return np.sqrt(x**2+y**2+z**2),np.arctan2(np.sqrt(x**2+y**2),z)
def totalN10State(x,y,z,n):
    return 1/np.sqrt(2*np.pi)*azimuthalN10State(np.arctan2(np.sqrt(x**2+y**2),z))*radialN10State(np.sqrt(x**2+y**2+z**2),n)
def total100State(x,y,z):
    return 1/(np.sqrt(np.pi)*a_0**(3/2))*np.exp(-np.sqrt(x**2+y**2+z**2)/a_0)
def integrand(r,n):
    return radialN10State(r,n)*radial100State(r)*r**3
def integrateRadialComponent(n,upper_bound,step_size):
    return np.sum(integrand(np.linspace(step_size//2,upper_bound+step_size/2,int(upper_bound//step_size)+1),n),axis=0)*step_size
def prefactor(k): # includes azimuthal 1/sqrt(3) integral term, 1/(E_n-E_k) term = 1/(e_1*(1-1/k**2)), and the radial integral R_n1^* R_10 r^3 dr
    # return 1
    return 1/(np.sqrt(3)*(e_1*(1-1/k**2)))*integrateRadialComponent(k,500*a_0,a_0/25)
def totalSumK10State(x,y,z,k):
    # return [(prefactor(n)*totalN10State(x,y,z,n))**2 for n in range(2,k+1)]
    # [print(n,prefactor(n))for n in range(2,k+1)]
    return np.sum([(prefactor(n)*totalN10State(x,y,z,n))**2 for n in range(2,k+1)],axis=0)
cubeLeng = 4 # this is how many radii out to go (a_0)
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
alpha_value = 0.07
s_value = 12
# ax = fig.add_subplot(122, projection='3d')
# colo = (totalN10State(xx,yy,zz,3))**2
# # colo = (total100State(xx,yy,zz))**2
# colo1=np.max(colo)
# colo = colo/np.max(colo)
# color_map = cm.ScalarMappable(cmap=cm.hsv)
# color_map.set_array(colo)
# img = ax.scatter(xx, yy, zz, marker='s',s=s_value,c=cm.hsv(colo),alpha=alpha_value)
# # plt.colorbar(color_map)
# # ax.set_title("3D Heatmap")
# ax.set_xticks(np.linspace(-cubeLeng*a_0,cubeLeng*a_0,5))
# ax.set_xticklabels(map(str,range(-cubeLeng,cubeLeng+1,cubeLeng//2)))
# ax.set_yticks(np.linspace(-cubeLeng*a_0,cubeLeng*a_0,5))
# ax.set_yticklabels(map(str,range(-cubeLeng,cubeLeng+1,cubeLeng//2)))
# ax.set_zticks(np.linspace(-cubeLeng*a_0,cubeLeng*a_0,5))
# ax.set_zticklabels(map(str,range(-cubeLeng,cubeLeng+1,cubeLeng//2)))
# ax.set_xlabel('X (a_0)')
# ax.set_ylabel('Y (a_0)')
# ax.set_zlabel('Z (a_0)')
# print(max(total100State(xx,yy,zz)**2))
# axe_field_mag = plt.axes([0.25, 0.9, 0.65, 0.03])
# e_field_mag = Slider(axe_field_mag,'Field Strength',0,1,0.5)
ax = fig.add_subplot(111, projection='3d')
# print(max((totalSumK10State(xx,yy,zz,12))**1))
# colo = total100State(xx,yy,zz)**2+e*e_field_mag.val*(totalSumK10State(xx,yy,zz,9))**1
colo = total100State(xx,yy,zz)**2+0*(totalSumK10State(xx,yy,zz,9))**1
# colo = (totalN10State(xx,yy,zz,9))**2
colo = colo/np.max(colo)
color_map = cm.ScalarMappable(cmap=cm.hsv)
color_map.set_array(colo)
img = ax.scatter(xx, yy, zz, marker='s',s=s_value,c=cm.hsv(colo),alpha=alpha_value)
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
def update(val):
    colo = (1-val)*total100State(xx,yy,zz)**2/max(total100State(xx,yy,zz)**2)+val*(totalSumK10State(xx,yy,zz,12))/max(totalSumK10State(xx,yy,zz,12))
    # colo = (1-e_field_mag.val)*total100State(xx,yy,zz)**2+e_field_mag.val*totalN10State(xx,yy,zz,12)**2
    colo = colo/np.max(colo)
    img.set_color(cm.hsv(colo))
# e_field_mag.on_changed(update)
# plt.colorbar(color_map)
animation = matplotlib.animation.FuncAnimation(fig,func=update,frames=np.linspace(0,0.25,26),interval=5,repeat=False)
animation.save('animation1.gif')
# displaying plot
plt.show()