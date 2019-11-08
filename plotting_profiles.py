import numpy as np
from pylab import *

out = np.loadtxt('multipoles.out')

r = out[:,0]
gt0 = out[:,1]
gt2 = out[:,2]
gx2 = out[:,3]

gt0_off = out[:,4]
gt_off  = out[:,5]
gt_cos_off = out[:,6]
gx_sin_off = out[:,7]

pcc = 0.8

gt = pcc*gt0 + (1-pcc)*gt0_off
gt_e = pcc*gt0 + (1-pcc)*gt_off

plt.figure()
plt.title('Gamma_t')
plt.plot(r,gt0,label = 'Only monopole')
plt.plot(r,gt,label = 'Monopole + misscentred (e=0)')
plt.plot(r,gt_e,label = 'Monopole + misscentred (e=0.25)')
plt.xscale('log')
plt.ylabel(r'$\Gamma_t0$')
plt.xlabel('r [Mpc]')
plt.legend()
plt.show()

# FOR ellip = 0.25

gt_cos = pcc*gt2 + (1-pcc)*gt_cos_off

plt.figure()
plt.title('Gamma_t_cos(2theta)')
plt.plot(r,gt2,label = 'Only quadrupole')
plt.plot(r,gt_cos,label = 'quadrupole + misscentred')
plt.xscale('log')
plt.ylabel(r'$\Gamma_t2$')
plt.xlabel('r [Mpc]')
plt.legend()
plt.show()

gx_sin = pcc*gx2 + (1-pcc)*gx_sin_off

plt.figure()
plt.title('Gamma_x_cos(2theta)')
plt.plot(r,gx2,label = 'Only quadrupole')
plt.plot(r,gx_sin,label = 'quadrupole + misscentred')
plt.xscale('log')
plt.ylabel(r'$\Gamma_x2$')
plt.xlabel('r [Mpc]')
plt.legend()
plt.show()


print 'END OF PROGRAM'
