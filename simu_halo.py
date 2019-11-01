import numpy as np
from pylab import *
from profiles_fit import *
from astropy.cosmology import LambdaCDM
from make_profile import *

cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

#parameters

cvel = 299792458;   # Speed of light (m.s-1)
G    = 6.670e-11;   # Gravitational constant (m3.kg-1.s-2)
pc   = 3.085678e16; # 1 pc (m)
Msun = 1.989e30 # Solar mass (kg)

############

rx = np.arange(-1.,1.,0.01)
ry = np.arange(-1.,1.,0.01)

x,y = np.meshgrid(rx,ry)

r = np.sqrt(x**2 + y**2)

mask = (r>0.1)

x,y = x[mask],y[mask]
r   = r[mask]


f = 1.

ab = np.sqrt((x**2 + (f*y)**2)/f)

M200 = 1.e14
z    = 0.2
zs   = 0.35

H        = cosmo.H(z).value/(1.0e3*pc) #H at z_pair s-1 
roc      = (3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_pair (kg.m-3)
roc_mpc  = roc*((pc*1.0e6)**3.0)
h        = cosmo.H(0).value/100.

disp = disp_sis(M200,H)
r200 = r200_nfw(M200,roc_mpc)

Ssis = sis_profile_sigma(ab,disp)

# plt.figure()
# plt.title('SIS')
# plt.scatter(x,y,c=np.log10(Ssis),vmin=3.0,vmax=2.0)
# plt.colorbar()

Snfw = NFW_profile_sigma((ab,roc_mpc,z,h),r200)
# plt.figure()
# plt.title('NFW')
# plt.scatter(x,y,c=np.log10(Snfw),vmin=3.0,vmax=2.0)
# plt.colorbar()

Dl    = cosmo.angular_diameter_distance(z).value*1.e6*pc
Ds    = cosmo.angular_diameter_distance(zs).value*1.e6*pc
Dls   = cosmo.angular_diameter_distance_z1z2(z,zs).value*1.e6*pc

Sc = ((((cvel**2.0)/(4.0*np.pi*G*Dl))*(1./(Dls/Ds)))*(pc**2/Msun))
sigma_c = np.zeros(len(r))
sigma_c.fill(Sc)

ex      = np.random.normal(0., 0.28, len(r))

profile_sis = shear_profile_log(100.,1000.,r*1.e3,Ssis/Sc,ex,np.ones(len(r)),np.zeros(len(r)),sigma_c,20)

R=profile_sis[0]/1.0e3 #r en Mpc
	
shear_sis  = profile_sis[1]
cero   = profile_sis[2]

profile_nfw = shear_profile_log(100.,1000.,r*1.e3,Snfw/Sc,ex,np.ones(len(r)),np.zeros(len(r)),sigma_c,20)

R=profile_sis[0]/1.0e3 #r en Mpc
	
shear_nfw  = profile_nfw[1]
cero   = profile_sis[2]



