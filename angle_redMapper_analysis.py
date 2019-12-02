import numpy as np
# from pylab import *
from astropy.io import fits
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70., Om0=0.3, Ode0=0.7)

folder = '/home/eli/Documentos/Astronomia/posdoc/halo-elongation/redMapper/'

members  = fits.open(folder+'redmapper_dr8_public_v6.3_members.fits')
clusters = fits.open(folder+'redmapper_dr8_public_v6.3_catalog.fits')
angles   = fits.open(folder+'angles_redMapper.fits')[1].data

ides  = members[1].data['ID']
R_cen = members[1].data['R']#[mask]
ra    = members[1].data['RA']#[mask]
ra[ra > 275] = ra[ra>275] - 360.
dec   = members[1].data['DEC']#[mask]
ID,c = np.unique(ides,return_counts=True)

ID_c  = clusters[1].data['ID']
zspec = clusters[1].data['Z_SPEC']
zlambda = clusters[1].data['Z_LAMBDA']
zc = zspec
zc[zc<0] = zlambda[zc<0]
Lambda = clusters[1].data['LAMBDA']


D_ang    = np.array(cosmo.angular_diameter_distance(zc))
kpcscale = D_ang*(((1.0/3600.0)*np.pi)/180.0)*1000.0

mcen = R_cen == 0.
RA0  = ra[mcen]
DEC0 = dec[mcen]

t     = angles['theta']
twl   = angles['theta_wlum']
tpwdl = angles['theta_pcut_wdl']
tp    = angles['theta_pcut']
tpwl  = angles['theta_pcut_wlum']
tpwd  = angles['theta_pcut_wd']

e     = angles['e']
ewl   = angles['e_wlum']
epwdl = angles['e_pcut_wdl']
ep    = angles['e_pcut']
epwl  = angles['e_pcut_wlum']
epwd  = angles['e_pcut_wd']

