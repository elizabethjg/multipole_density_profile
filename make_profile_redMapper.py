import numpy as np
# from pylab import *
from astropy.io import fits
#parameters
from astropy.cosmology import LambdaCDM
from maria_func import *
cosmo = LambdaCDM(H0=70., Om0=0.3, Ode0=0.7)

folder = '/mnt/clemente/lensing/redMaPPer/'

cfht = fits.open(folder+'gx_CFHT_redMapper.fits')[1].data
kids = fits.open(folder+'gx_KiDS_redMapper.fits')[1].data
cs82 = fits.open(folder+'gx_CS82_redMapper.fits')[1].data

print '---- EXTRACT DATA -------'

ID = np.concatenate((cfht.ID,kids.ID,cs82.ID))

ra     = np.concatenate((cfht.RAJ2000,kids.RAJ2000,cs82.RAJ2000))
dec    = np.concatenate((cfht.DECJ2000,kids.DECJ2000,cs82.DECJ2000))
	
e1     = np.concatenate((cfht.e1,kids.e1,cs82.e1))
e2     = np.concatenate((cfht.e2,kids.e2,cs82.e2))
	
zlambda = np.concatenate((cfht.Z_LAMBDA,kids.Z_LAMBDA,cs82.Z_LAMBDA))
zspec   = np.concatenate((cfht.Z_SPEC,kids.Z_SPEC,cs82.Z_SPEC))
Zl      = zspec
Zl[Zl<0] = zlambda[Zl<0]
	
	
ALFA0  = np.concatenate((cfht.RA,kids.RA,cs82.RA))
DELTA0 = np.concatenate((cfht.DEC,kids.DEC,cs82.DEC))

dls  = np.concatenate((cfht.DLS,kids.DLS,cs82.DLS))
ds   = np.concatenate((cfht.DS,kids.DS,cs82.DS))
dl   = np.concatenate((cfht.DL,kids.DL,cs82.DL))

peso = np.concatenate((cfht.weight,kids.weight,cs82.weight))
m    = np.concatenate((cfht.m,kids.m,cs82.m))


# ides,index = np.unique(ID,return_index=True)

# NCUM = len(ides)

# print 'cantidad de lentes',NCUM


KPCSCALE   = dl*(((1.0/3600.0)*np.pi)/180.0)*1000.0
BETA_array = dls/ds
beta       = BETA_array.mean()

Dl = dl*1.e6*pc
sigma_c = (((cvel**2.0)/(4.0*np.pi*G*Dl))*(1./BETA_array))*(pc**2/Msun)


print 'BETA.mean',beta

SIGMAC = (((cvel**2.0)/(4.0*np.pi*G*Dl.mean())))*(pc**2/Msun)

print 'SIGMA_C', SIGMAC


rads, theta, test1,test2 = eq2p2(np.deg2rad(ra),
					np.deg2rad(dec),
					np.deg2rad(ALFA0),
					np.deg2rad(DELTA0))


#Correct polar angle for e1, e2
theta = theta+np.pi/2.

#get tangential ellipticities 
et = (e1*np.cos(2*theta)+e2*np.sin(2*theta))*sigma_c
#get cross ellipticities
ex = (-e1*np.sin(2*theta)+e2*np.cos(2*theta))*sigma_c


r=np.rad2deg(rads)*3600*KPCSCALE


zmean    = (Zl).mean()
zdisp    = (Zl).std()
H        = cosmo.H(zmean).value/(1.0e3*pc) #H at z_pair s-1 
roc      = (3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_pair (kg.m-3)
roc_mpc  = roc*((pc*1.0e6)**3.0)
D_ang    = cosmo.angular_diameter_distance(zmean)
kpcscale = D_ang*(((1.0/3600.0)*np.pi)/180.0)*1000.0


print '---------------------------------------------------------'
print '             COMPUTING THE SHEAR PROFILES                '
print '========================================================='

profile = shear_profile_log(RIN,ROUT,r,et,ex,peso,m,sigma_c,ndots,
                            booterror_flag=False)



R=profile[0]/1.0e3 #r en Mpc

shear  = profile[1]
cero   = profile[2]
BIN    = profile[4]
err_et = profile[5]
err_ex = profile[6]

print 'Now is trying to fit a NFW profile...'

try:
	nfw          = NFW_stack_fit(R[:BIN],shear[:BIN],err_et[:BIN],zmean,roc)
except:
	nfw          = [-999.,0.,-100.,[0.,0.],[0.,0.],-999.,0.]

c          = nfw[5]
CHI_nfw    = nfw[2]
R200       = nfw[0]
error_R200 = nfw[1]
x2         = nfw[3]
y2         = nfw[4]
