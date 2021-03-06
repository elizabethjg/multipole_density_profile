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

r  = np.logspace(np.log10(0.05),np.log10(5.),30)

'''
out = multipole_shear_parallel(r,M200=4.14e14,z=0.21,
                              ellip=0.,misscentred=True,
                              ncores=20)

save = np.array([r,out['Gt0'],out['Gt2'],out['Gx2'],out['Gt0_off'],
               out['Gt_off'],out['Gt_off_cos'],out['Gx_off_sin']]).T
# save = np.array([r,out['Gt0'],out['Gt2'],out['Gx2']]).T

np.savetxt('../multipoles_bin4.out',save)


'''

for q in np.arange(0.1,1.,0.2):

     print '############################'
     print 'COMPUTING FOR q ',q
     print '----------------------------'
     
     ellip = (1.- q)/(1. + q)
     
     out = multipole_shear_parallel(r,M200=M,z=0.2,
                                   ellip=ellip,misscentred=True,
                                   verbose=False,ncores=30)
     
     save = np.array([r,out['Gt0'],out['Gt2'],out['Gx2'],out['Gt0_off'],
                    out['Gt_off'],out['Gt_off_cos'],out['Gx_off_sin']]).T
     # save = np.array([r,out['Gt0'],out['Gt2'],out['Gx2']]).T
     
     np.savetxt('../multipoles_'+str('%2.1f' % q)+'.out',save)


ellip = 0.25

for logM in np.arange(12.5,15.5,0.5):
     
     print '############################'
     print 'COMPUTING FOR MASS ',logM
     print '----------------------------'
     
     M = 10**logM
     
     out = multipole_shear_parallel(r,M200=M,z=0.2,
                                   ellip=ellip,misscentred=True,
                                   verbose=False,ncores=30)
     
     save = np.array([r,out['Gt0'],out['Gt2'],out['Gx2'],out['Gt0_off'],
                    out['Gt_off'],out['Gt_off_cos'],out['Gx_off_sin']]).T
 
     
     np.savetxt('../multipoles_'+str(logM)+'.out',save)

#'''

print 'END OF PROGRAM'
