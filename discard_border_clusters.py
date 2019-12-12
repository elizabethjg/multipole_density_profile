import numpy as np
from astropy.io import fits
from astropy.cosmology import LambdaCDM
from ang_sep import *
cosmo = LambdaCDM(H0=70., Om0=0.3, Ode0=0.7)

folder = '/home/eli/Documentos/Astronomia/posdoc/halo-elongation/redMapper/'

members  = fits.open(folder+'redmapper_dr8_public_v6.3_members.fits')[1].data
clusters = fits.open(folder+'redmapper_dr8_public_v6.3_catalog.fits')[1].data

ides = members.ID
ID,c = np.unique(ides,return_counts=True)

R_cen = members.R
ra    = members.RA
ra[ra > 275] = ra[ra>275] - 360.
dec   = members.DEC

bines_ra  = np.linspace(ra.min(),ra.max(),num=200)
bines_dec = np.linspace(dec.min(),dec.max(),num=200)

zspec = clusters.Z_SPEC
zlambda = clusters.Z_LAMBDA
zc = zspec
zc[zc<0] = zlambda[zc<0]
D_ang    = np.array(cosmo.angular_diameter_distance(zc))
kpcscale = D_ang*(((1.0/3600.0)*np.pi)/180.0)*1000.0
KPCSCALE = np.repeat(kpcscale,c)


ID_border = np.array([])


# plt.plot(ra[mask],dec[mask],'k,')

for BIN in np.arange(len(bines_ra)-1):

     ramin = bines_ra[BIN]
     ramax = bines_ra[BIN+1]
     mask = (ra >= ramin)*(ra<ramax)
     
     if mask.sum() != 0.:
          
          # plt.plot(ra[mask],dec[mask])
     
          j = np.argmin(dec[mask])
          minra  = ra[mask][j]
          mindec = dec[mask][j]
          
          j = np.argmax(dec[mask])
          maxra  = ra[mask][j]
          maxdec = dec[mask][j]
          
          
          mcen = R_cen[mask] == 0.
          RA0  = ra[mask][mcen]
          DEC0 = dec[mask][mcen]
          
          distmin = ang_sep(minra,mindec,RA0,DEC0)*3600*(KPCSCALE[mask][mcen])
          distmax = ang_sep(maxra,maxdec,RA0,DEC0)*3600*(KPCSCALE[mask][mcen])
          
          mdist = (distmin < 2000.)+(distmax < 2000.)
          
          ID_border = np.append(ID_border,ides[mask][mcen][mdist])

     decmin = bines_dec[BIN]
     decmax = bines_dec[BIN+1]
     mask = (dec >= decmin)*(dec<decmax)


     if mask.sum() != 0.:
     
          j = np.argmin(ra[mask])
          minra  = ra[mask][j]
          mindec = dec[mask][j]
          
          j = np.argmax(ra[mask])
          maxra  = ra[mask][j]
          maxdec = dec[mask][j]
          
          j = np.argmin(abs(ra-60.)[mask])
          medra1  = ra[mask][j]
          meddec1 = dec[mask][j]

          j = np.argmin(abs(ra-100.)[mask])
          medra2  = ra[mask][j]
          meddec2 = dec[mask][j]
          
          mcen = R_cen[mask] == 0.
          RA0  = ra[mask][mcen]
          DEC0 = dec[mask][mcen]
          
          distmin  = ang_sep(minra,mindec,RA0,DEC0)*3600*(KPCSCALE[mask][mcen])
          distmax  = ang_sep(maxra,maxdec,RA0,DEC0)*3600*(KPCSCALE[mask][mcen])
          distmed1 = ang_sep(medra1,meddec1,RA0,DEC0)*3600*(KPCSCALE[mask][mcen])
          distmed2 = ang_sep(medra2,meddec2,RA0,DEC0)*3600*(KPCSCALE[mask][mcen])
          
          mdist = (distmin < 2000.)+(distmax < 2000.)+(distmed1 < 2000.)+(distmed2 < 2000.)
          
          ID_border = np.append(ID_border,ides[mask][mcen][mdist])

ID_border = np.unique(ID_border)

border = np.in1d(ides,ID_border)
plt.plot(ra,dec,'k,')
plt.plot(ra[border],dec[border],'.')

np.savetxt(folder+'redMapperID_border.list',ID_border,fmt=['%7i'])
