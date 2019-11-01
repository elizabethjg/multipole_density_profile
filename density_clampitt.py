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

# a  = np.logspace(np.log10(0.01),np.log10(5.),10)
a  = np.linspace(0.01,5.,10)
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

q = 0.6
ellip = (1.- q**2)/(1. + q**2) 

print '---------------------'
print '      MONOPOLE       '
print '---------------------'

P = multipole_clampitt(r,theta,f=q)

DS0 = P['monopole']

DS0_1 = -1*DS0*np.cos(2*theta)
DS0_2 = -1*DS0*np.sin(2*theta)

DS0x1 = np.zeros(len(DS0_1))
DS0y1 = np.zeros(len(DS0_1))
DS0x2 = np.zeros(len(DS0_1))
DS0y2 = np.zeros(len(DS0_1))

DS0x1[DS0_1 > 0.] = (DS0_1[DS0_1 > 0.])*np.sign(np.sin(theta[DS0_1 > 0.]))
DS0y1[DS0_1 < 0.] = (DS0_1[DS0_1 < 0.])*np.sign(np.cos(theta[DS0_1 < 0.]))

DS0x2 = DS0_2*np.cos(np.pi/4.)*-1*np.sign(np.cos(theta))
DS0y2 = np.abs(DS0_2*np.sin(np.pi/4.))*-1*np.sign(np.cos(theta))


DS0x = DS0x1 + DS0x2
DS0y = DS0y1 + DS0y2

DS0norm  = np.sqrt(DS0x**2 + DS0y**2)
DS0norm1 = np.sqrt(DS0x1**2 + DS0y1**2)
DS0norm2 = np.sqrt(DS0x2**2 + DS0y2**2)

DS0x_normed = DS0x / DS0norm
DS0y_normed = DS0y / DS0norm

f, ax = plt.subplots(1, 3, sharey=True,figsize=(14,3))

ax[0].quiver(x,y,DS0x_normed,DS0y_normed,DS0norm)
ax[0].set_xlabel('x [Mpc]')
ax[0].set_ylabel('y [Mpc]')
ax[1].quiver(x,y,DS0x2,DS0y2,DS0norm2)
ax[1].set_xlabel('x [Mpc]')
ax[1].set_title('Monopole')
ax[2].quiver(x,y,DS0x1,DS0y1,DS0norm1)
ax[2].set_xlabel('x [Mpc]')

print '---------------------'
print '      QUADRUPOLE     '
print '---------------------'


DS_1 = (ellip/4.)*((2.*P['I1'] - P['quadrupole'])*np.cos(4*theta) +
                   2.*P['I2'] - P['quadrupole'])

DS_2 = (ellip/4.)*(2.*P['I1'] - P['quadrupole'])*np.sin(4*theta)


DS1x = np.zeros(len(DS_1))
DS1y = np.zeros(len(DS_1))
DS2x = np.zeros(len(DS_1))
DS2y = np.zeros(len(DS_1))

DS1x[DS_1>0]  = (DS_1[DS_1>0])#*np.sign(np.sin(theta[DS_1 > 0.]))
DS1y[DS_1<0]  = np.abs(DS_1[DS_1<0])#*-1*np.sign(np.cos(theta[DS_1 < 0.]))

DS2x  = DS_2*np.cos(np.pi/4.)#*-1*np.sign(np.cos(theta))
DS2y  = np.abs(DS_2)*np.sin(np.pi/4.)#*-1*np.sign(np.sin(theta))

DSx = DS1x + DS2x
DSy = DS1y + DS2y

DSnorm  = np.sqrt(DSx**2 + DSy**2)
DSnorm1 = np.sqrt(DS1x**2 + DS1y**2)
DSnorm2 = np.sqrt(DS2x**2 + DS2y**2)

DSx_normed = DSx / DSnorm
DSy_normed = DSy / DSnorm

# plt.figure()
# plt.quiver(x,y,DSx,DSy)

f, ax = plt.subplots(1, 3, sharey=True,figsize=(14,3))
ax[0].quiver(x,y,DSx_normed,DSy_normed,DSnorm)
ax[0].set_xlabel('x [Mpc]')
ax[0].set_ylabel('y [Mpc]')
ax[1].set_title('Quadrupole')
ax[1].quiver(x,y,DS2x,DS2y,DSnorm2)
ax[1].set_xlabel('x [Mpc]')
ax[2].quiver(x,y,DS1x,DS1y,DSnorm1)
ax[2].set_xlabel('x [Mpc]')


print '----------------------------'
print '     QUADRUPOLE PROFILE     '
print '----------------------------'
# Adhikari et al. (2015) - eq 14 y 15

DSt = (P['I1'] + P['I2'] - P['quadrupole'])*(ellip/2.)
DSx = (P['I1'] - P['I2'])*(ellip/2.)

plt.figure()
plt.plot(r,DSt,label=r'$\Delta \Sigma_+$')
plt.plot(r,DSx,label=r'$\Delta \Sigma_x$')
plt.xscale('log')
plt.axis([0.1,5,-400,400])
plt.legend()


# '''
'''
print '----------------------------'
print '   QUADRUPOLE MISCENTRED    '
print '----------------------------'


DS_1_off = (ellip/4.)*((2.*I1r_off - quadrupole_off_r)*np.cos(4*theta) +
               2.*I2r_off - quadrupole_off_r)

DS_2_off = (ellip/4.)*(2.*I1r_off - quadrupole_off_r)*np.sin(4*theta)


DS1x_off = np.zeros(len(DS_1_off))
DS1y_off = np.zeros(len(DS_1_off))
DS2x_off = np.zeros(len(DS_1_off))
DS2y_off = np.zeros(len(DS_1_off))

DS1x_off[DS_1_off>0]  = (DS_1_off[DS_1_off>0])#*np.sign(np.sin(theta[DS_1 > 0.]))
DS1y_off[DS_1_off<0]  = np.abs(DS_1_off[DS_1_off<0])#*-1*np.sign(np.cos(theta[DS_1 < 0.]))

DS2x_off  = DS_2_off*np.cos(np.pi/4.)#*-1*np.sign(np.cos(theta))
DS2y_off  = np.abs(DS_2_off)*np.sin(np.pi/4.)#*-1*np.sign(np.sin(theta))

DSx_off = DS1x_off + DS2x_off
DSy_off = DS1y_off + DS2y_off

DSnorm_off  = np.sqrt(DSx_off**2 + DSy_off**2)
DSnorm1_off = np.sqrt(DS1x_off**2 + DS1y_off**2)
DSnorm2_off = np.sqrt(DS2x_off**2 + DS2y_off**2)

DSx_normed_off = DSx_off / DSnorm_off
DSy_normed_off = DSy_off / DSnorm_off

# plt.figure()
# plt.quiver(x,y,DSx,DSy)

f, ax = plt.subplots(1, 3, sharey=True,figsize=(14,3))
# ax[0].quiver(x,y,DSx_normed_off,DSy_normed_off,np.log10(DSnorm_off))
ax[0].quiver(x,y,DSx_off,DSy_off,DSnorm_off)
ax[1].quiver(x,y,DS2x_off,DS2y_off,DSnorm2_off)
ax[2].quiver(x,y,DS1x_off,DS1y_off,DSnorm1_off)


DSt_off = (I1r_off + I2r_off - quadrupole_off_r)*(e/2.)
DSx_off = (I1r_off - I2r_off)*(e/2.)


'''




'''
r = np.linspace(0.1,5,500)

I1r = vecI1(r)
I2r = vecI2(r)

DS_4t = (ellip/4.)*(2.*I1r - monopole(r)*etha0(r))
DS_c  = (ellip/4.)*(2.*I2r - monopole(r)*etha0(r))

plt.figure()
plt.plot([0,10],[0,0],'k--')
plt.plot(r,DS_4t)
plt.plot(r,DS_c)
plt.axis([0.1,5.,-45,10])
plt.show()

'''
