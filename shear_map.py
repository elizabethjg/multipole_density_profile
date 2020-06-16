import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
import numpy as np
from matplotlib import *
from multipoles_shear import *
from astropy.io import fits
import corner
from profiles_fit import *
import os
#parameters

M200 = 2.e14
e    = 0.5
z    = 0.25

############  MAKING A GRID

# a  = np.logspace(np.log10(0.01),np.log10(5.),10)
a  = np.arange(-1,1.3,0.25)
# a  = np.append(a,-1.*a)

x,y = np.meshgrid(a,a)

x = x.flatten()
y = y.flatten()

r = np.sqrt(x**2 + y**2)

mask = (r>0.2)

x,y = x[mask],y[mask]
r   = r[mask]
theta  = np.arctan2(y,x)

j   = argsort(r)
r   = r[j]
theta  = theta[j]
x,y = x[j],y[j]

# COMPUTE COMPONENTS

out = multipole_shear_parallel(r,M200=M200,z=z,ellip=e)


f, ax = plt.subplots(1,3, figsize=(12,4), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)
######################
# MONOPOLE
######################

g0_1 = -1.*out['Gt0']*np.cos(2.*theta)
g0_2 = -1.*out['Gt0']*np.sin(2.*theta)

g0_x = np.zeros(len(r))
g0_y = np.zeros(len(r))

g0_x[g0_1 > 0.] -= g0_1[g0_1 > 0.]*np.sign(np.sin(theta[g0_1 > 0.])*np.cos(theta[g0_1 > 0.])) 
g0_x += g0_2*np.cos(np.pi/4.)

g0_y[g0_1 < 0.] -= g0_1[g0_1 < 0.]
g0_y += np.abs(g0_2*np.sin(np.pi/4.))

g0_norm = np.sqrt(g0_x**2 + g0_y**2)
mq = ax[0].quiver(x,y,g0_x/g0_norm,g0_y/g0_norm,g0_norm,pivot='mid',headlength=0, headwidth = 1,norm=matplotlib.colors.LogNorm())

######################
# QUADRUPOLE
######################

gt2 = out['Gt2']*np.cos(2.*theta)
gx2 = out['Gx2']*np.sin(2.*theta)

g2_1 = -1.*gt2*np.cos(2.*theta) + gx2*np.sin(2.*theta)
g2_2 = -1.*gt2*np.sin(2.*theta) - gx2*np.cos(2.*theta)

g2_x = np.zeros(len(r))
g2_y = np.zeros(len(r))

g2_x[g2_1 > 0.] -= g2_1[g2_1 > 0.]*np.sign(np.sin(theta[g2_1 > 0.])*np.cos(theta[g2_1 > 0.])) 
g2_x += g2_2*np.cos(np.pi/4.)

g2_y[g2_1 < 0.] -= g2_1[g2_1 < 0.]
g2_y += np.abs(g2_2*np.sin(np.pi/4.))

g2_norm = np.sqrt(g2_x**2 +g2_y**2)
qq = ax[1].quiver(x,y,g2_x/g2_norm,g2_y/g2_norm,g2_norm,pivot='mid',headlength=0, headwidth = 1,norm=matplotlib.colors.LogNorm())

######################
# TOTAL SHEAR MAP
######################
gt = out['Gt0'] + out['Gt2']*np.cos(2.*theta)
gx = out['Gx2']*np.sin(2.*theta)

g_1 = -1.*gt*np.cos(2.*theta) + gx*np.sin(2.*theta)
g_2 = -1.*gt*np.sin(2.*theta) - gx*np.cos(2.*theta)

g_x = np.zeros(len(r))
g_y = np.zeros(len(r))

g_x[g_1 > 0.] -= g_1[g_1 > 0.]*np.sign(np.sin(theta[g_1 > 0.])*np.cos(theta[g_1 > 0.])) 
g_x += g_2*np.cos(np.pi/4.)

g_y[g_1 < 0.] -= g_1[g_1 < 0.]
g_y += np.abs(g_2*np.sin(np.pi/4.))

g_norm = np.sqrt(g_x**2 +g_y**2)

tq = ax[2].quiver(x,y,g_x/g_norm,g_y/g_norm,g_norm,pivot='mid',headlength=0, headwidth = 1,norm=matplotlib.colors.LogNorm())

plt.axis([-1.2,1.2,-1.2,1.2])

mq.set_clim(20,80)
qq.set_clim(20,80)
tq.set_clim(20,80)
# cbar = f.colorbar(ticks=[10,20,30,40])
# cbar.ax.set_yticklabels(['10','20','30'])
