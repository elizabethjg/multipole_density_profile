import numpy as np
from pylab import *
from astropy.cosmology import LambdaCDM
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

r  = np.logspace(np.log10(0.05),np.log10(1.5),4)

ellip = (1.- q)/(1. + q)
out = multipole_shear_parallel(r,M200=M,z=0.2,zs=0.6,ellip=ellip,misscentred=True,ncores=2)

save = np.array([r,out['Gt0'],out['Gt2'],out['Gx2'],out['Gt0_off'],out['Gt_off'],out['Gx_off']]).T
# save = np.array([r,out['Gt0'],out['Gt2'],out['Gx2']]).T

np.savetxt('../multipoles.out',save)

print 'END OF PROGRAM'
