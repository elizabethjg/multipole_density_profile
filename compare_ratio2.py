import numpy as np
from pylab import *
from profiles_fit import *
from astropy.cosmology import LambdaCDM
from make_profile import *
from scipy.misc import derivative
from scipy import integrate
from multipoles_shear import *

#parameters

cvel = 299792458;   # Speed of light (m.s-1)
G    = 6.670e-11;   # Gravitational constant (m3.kg-1.s-2)
pc   = 3.085678e16; # 1 pc (m)
Msun = 1.989e30 # Solar mass (kg)

############  MAKING A GRID

logM = 14.
M    = 10**logM
q    = 0.6

a  = np.logspace(np.log10(0.05),np.log10(1.5),15)
a  = np.append(a,-1.*a)

x,y = np.meshgrid(a,a)

x = x.flatten()
y = y.flatten()

r = np.sqrt(x**2 + y**2)

mask = (r>0.01)

x,y = x[mask],y[mask]
r   = r[mask]
theta  = np.arctan2(y,x)

j   = argsort(r)
r   = r[j]
theta  = theta[j]
x,y = x[j],y[j]

# q = 0.6

f, ax = plt.subplots(2, 2, figsize=(10,8))
f2, ax2 = plt.subplots(1, 2, figsize=(14,7))

print '----------------------------'
print '     QUADRUPOLE PROFILE     '
print '----------------------------'


# for q in np.arange(0.1,1.,0.2):
for logM in np.arange(12.5,15,0.5):
     
     print q
     print logM
     
     M = 10**logM
     
     Pa = multipole_clampitt(r,theta,M200=M,z=0.2,zs=0.6)
     ellip = (1.- q**2)/(2.*(1. + q**2))
     
     # Adhikari et al. (2015) - eq 14 y 15
     
     DSt = (Pa['I1'] + Pa['I2'] - Pa['quadrupole'])*(ellip)
     DSx = (Pa['I1'] - Pa['I2'])*(ellip)
     
     ax[0,0].plot([0,5],[0,0],'k--')
     ax[0,0].plot(r,DSt,str(q/2.))
     ax[0,0].set_xscale('log')
     ax[0,0].set_ylabel(r'$\Delta \Sigma^c_+$')
     ax[0,1].plot([0,1.5],[0,0],'k--')
     ax[0,1].plot(r,DSx,str(q/2.))
     ax[0,1].set_ylabel(r'$\Delta \Sigma^c_x$')
     ax[0,1].set_xscale('log')
     ax[0,1].set_xlabel('r [Mpc]')
     
     ax2[0].plot([0,1.5],[0,0],'k--')
     ax2[1].plot([0,1.5],[0,0],'k--')
     ax2[0].plot(r,DSt,'C0',alpha=1.-(q/2.))
     ax2[0].set_ylabel(r'$\Delta \Sigma_+$')
     ax2[1].plot(r,DSx,'C0',alpha=1.-(q/2.))
     ax2[1].set_ylabel(r'$\Delta \Sigma_x$')
     
     Pv = multipole_vanUitert(r,theta,M200=M,z=0.2,zs=0.6)
     # ellip = (1.- q)/(1. + q) 
     
     DS_t = ellip*((-6*Pv['psi2']/r**2) - 2.*Pv['monopole'] + Pv['quadrupole'])
     DS_x = ellip*((-6*Pv['psi2']/r**2) - 4.*Pv['monopole'])

     ax[1,0].plot([0,1.5],[0,0],'k--')
     ax[1,0].plot(r,DS_t,str(q/2.))
     ax[1,0].set_xscale('log')
     ax[1,0].set_ylabel(r'$\Delta \Sigma^{vU}_+$')
     ax[1,1].plot([0,1.5],[0,0],'k--')
     ax[1,1].plot(r,DS_x,str(q/2.))
     ax[1,1].set_ylabel(r'$\Delta \Sigma^{vU}_x$')
     ax[1,1].set_xscale('log')
     ax[1,1].set_ylabel('r [Mpc]')

     ax2[0].plot(r,DS_t,'C1',alpha=1.-(q/2.))
     ax2[1].plot(r,DS_x,'C1',alpha=1.-(q/2.))
     ax2[0].set_xscale('log')
     ax2[1].set_xscale('log')
     ax2[0].set_xlabel('r [Mpc]')
     ax2[1].set_xlabel('r [Mpc]')
     
plt.show()
