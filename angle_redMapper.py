import numpy as np
# from pylab import *
from multipoles_shear import *
from astropy.io import fits
#parameters

folder = '/home/eli/Documentos/Astronomia/posdoc/halo-elongation/redMapper/'

members  = fits.open(folder+'redmapper_dr8_public_v6.3_members.fits')
clusters = fits.open(folder+'redmapper_dr8_public_v6.3_catalog.fits')

ides  = members[1].data['ID']
RA    = members[1].data['RA']
DEC   = members[1].data['DEC']
R_cen = members[1].data['R']
Lum_r = 10.**(-0.4*members[1].data['MODEL_MAG_R'])

ID_c  = clusters[1].data['ID']

e1 = np.array([])
e2 = np.array([])

for j in range(len(ID_c)):
     mid  = ides == ID_c[j]
     ra   = RA[mid]
     dec  = DEC[mid]
     lum  = Lum_r[mid]
     mcen = R_cen[mid] == 0
     ra0  = ra[mcen]
     dec0 = dec[mcen]
     Q11  = np.sum((ra[~mcen]-ra0)**2*lum[~mcen])/np.sum(lum[~mcen])
     Q22  = np.sum((dec[~mcen]-dec0)**2*lum[~mcen])/np.sum(lum[~mcen])
     Q12  = np.sum((ra[~mcen]-ra0)*(dec[~mcen]-dec0)*lum[~mcen])/np.sum(lum[~mcen])
     e1 = np.append(e1,(Q11-Q22)/(Q11+Q22))
     e2 = np.append(e1,(2.*Q12)/(Q11+Q22))
