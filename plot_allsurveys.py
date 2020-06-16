import numpy as np
from matplotlib import *
from multipoles_shear import *
from astropy.io import fits

folder = '/mnt/clemente/lensing/redMaPPer/'

cfht     = fits.open(folder+'gx_CFHT_redMapper.fits')[1].data
kids     = fits.open(folder+'gx_KiDS_redMapper.fits')[1].data
cs82     = fits.open(folder+'gx_CS82_redMapper.fits')[1].data
rcsl     = fits.open(folder+'gx_RCSL_redMapper.fits')[1].data
angles   = fits.open(folder+'angles_redMapper_forprofile.fits')[1].data
clusters = fits.open(folder+'redmapper_dr8_public_v6.3_catalog.fits')[1].data
borderid = np.loadtxt(folder+'redMapperID_border.list')

mkids = (~np.in1d(kids.ID,cfht.ID))*(kids.Z_B < 0.9)
mrcsl = (~np.in1d(rcsl.ID,cfht.ID))*(~np.in1d(rcsl.ID,cs82.ID))*(~np.in1d(rcsl.ID,kids.ID))*(rcsl.Z_B < 1.3)
mcs82 = (~np.in1d(cs82.ID,cfht.ID))*(cs82.Z_B < 1.3)
mcfht = (cfht.Z_B < 1.3)

kids = kids[mkids]
cs82 = cs82[mcs82]
rcsl = rcsl[mrcsl]
cfht = cfht[mcfht]

RA = cs82.RAJ2000
DEC = cs82.RAJ2000
RA[RA < 0] = RA[RA<0] + 360.

plt.plot(clusters.RA,clusters.DEC,'ko')
plt.plot(cfht.RAJ2000,cfht.DECJ2000,'C3,')
plt.plot(kids.RAJ2000,kids.DECJ2000,'C4,')
plt.plot(RA,DEC,'C5,',label='CS82')
plt.plot(rcsl.RAJ2000,rcsl.DECJ2000,'C6,')
plt.plot(360,360,'C3o',label='CFTHLens')
plt.plot(360,360,'C4o',label='KiDS-450')
plt.plot(360,360,'C5o',label='CS82')
plt.plot(360,360,'C6o',label='RCSLens')
plt.legend()
plt.xlabel('R.A. [deg]')
plt.ylabel('Dec. [deg]')
plt.axis([0,360,-2,80])
plt.savefig('redmapper_surveys.png',bbox_inches='tight')
plt.savefig('/home/elizabeth/redmapper_surveys.png',bbox_inches='tight')
