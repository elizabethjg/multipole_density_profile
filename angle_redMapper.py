import numpy as np
# from pylab import *
from astropy.io import fits
#parameters
from astropy.cosmology import LambdaCDM
from astropy.wcs import WCS
cosmo = LambdaCDM(H0=70., Om0=0.3, Ode0=0.7)

wcs = WCS(naxis=2)
wcs.wcs.crpix = [0., 0.]
wcs.wcs.cdelt = [1./3600., 1./3600.]
wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

def momentos(dx,dy,w):
     
     Q11  = np.sum((dx**2)*w)/np.sum(w)
     Q22  = np.sum((dy**2)*w)/np.sum(w)
     Q12  = np.sum((dx*dy*2)*w)/np.sum(w)
     E1 = (Q11-Q22)/(Q11+Q22)
     E2 = (2.*Q12)/(Q11+Q22)
     e = np.sqrt(E1**2 + E2**2)
     theta = np.arctan2(E2,E1)/2.
     return e,theta
     

folder = '/home/eli/Documentos/Astronomia/posdoc/halo-elongation/redMapper/'

members  = fits.open(folder+'redmapper_dr8_public_v6.3_members.fits')
clusters = fits.open(folder+'redmapper_dr8_public_v6.3_catalog.fits')

ides  = members[1].data['ID']
ID,c = np.unique(ides,return_counts=True)

R_cen = members[1].data['R']
RA    = members[1].data['RA']
RA[RA > 275] = RA[RA>275] - 360.
DEC   = members[1].data['DEC']

P     = members[1].data['P']
zspec = clusters[1].data['Z_SPEC']
zlambda = clusters[1].data['Z_LAMBDA']
zc = zspec
zc[zc<0] = zlambda[zc<0]
zc = np.repeat(zc,c)
D_lum    = np.array(cosmo.luminosity_distance(zc))*1.e6
MAG_abs=members[1].data['MODEL_MAG_R']+5.0-5.0*np.log10(D_lum)

Lum_r = 10.**(-0.4*MAG_abs)

ID_c  = clusters[1].data['ID']

N = np.array([])

e = np.array([])
theta = np.array([])

e_lum = np.array([])
theta_lum = np.array([])

e_pcutdl = np.array([])
theta_pcutdl = np.array([])

e_pcut = np.array([])
theta_pcut = np.array([])

e_pcutl = np.array([])
theta_pcutl = np.array([])

e_pcutd = np.array([])
theta_pcutd = np.array([])

e_control = np.array([])
theta_control = np.array([])


X = np.array([])
Y = np.array([])

plots = np.linspace(1,len(ID_c),10)

f,ax = plt.subplots()
f2,ax2 = plt.subplots()

for j in ind.astype(int):
     mid  = ides == j
# for j in range(len(ID_c)):
     print j
     # mid  = ides == ID_c[j]
     N = np.append(N,mid.sum())
     ra   = RA[mid]
     dec  = DEC[mid]
     mcen = R_cen[mid] == 0
     
     ra0  = ra[mcen]
     dec0 = dec[mcen]
     
     wcs.wcs.crval = [ra0[0],dec0[0]]
     dx, dy = wcs.wcs_world2pix(ra, dec, 0)

     X = np.append(X,dx)
     Y = np.append(Y,dy)


     dx = dx[~mcen]
     dy = dy[~mcen]

     # control
     ellip, ang = momentos(ra[~mcen]-ra0,dec[~mcen]-dec0,np.ones(len(dx)))
     e_control = np.append(e_control,ellip)
     theta_control = np.append(theta_control,ang)     
     
     # all members
     ellip, ang = momentos(dx,dy,np.ones(len(dx)))
     e = np.append(e,ellip)
     theta = np.append(theta,ang)
     
     # weighted by luminosity
     wl    = Lum_r[mid][~mcen]
     ellip, ang = momentos(dx,dy,wl)
     e_lum = np.append(e_lum,ellip)
     theta_lum = np.append(theta_lum,ang)
          
     # p_cut
     pcut = P[mid][~mcen] > 0.5
     ellip, ang = momentos(dx[pcut],dy[pcut],np.ones(len(dx[pcut])))
     e_pcut = np.append(e_pcut,ellip)
     theta_pcut = np.append(theta_pcut,ang)
     
     # p_cut weighted by lum
     ellip, ang = momentos(dx[pcut],dy[pcut],wl[pcut])
     e_pcutl = np.append(e_pcutl,ellip)
     theta_pcutl = np.append(theta_pcutl,ang)    
     
     # p_cut weighted by distance
     wd   = 1./(dx**2 + dy**2)
     ellip, ang = momentos(dx[pcut],dy[pcut],wd[pcut])
     e_pcutd = np.append(e_pcutd,ellip)
     theta_pcutd = np.append(theta_pcutd,ang)         
     
     # p_cut weighted by distance and lum
     wd   = Lum_r[mid][~mcen]/(dx**2 + dy**2)
     ellip, ang = momentos(dx[pcut],dy[pcut],wd[pcut])
     e_pcutdl = np.append(e_pcutdl,ellip)
     theta_pcutdl = np.append(theta_pcutdl,ang)
     
     # plt.plot(ra,dec,'r.')
     
     ax.plot(ra-ra0,dec-dec0,'C0.')
     ax2.plot(dx,dy,'C1.')
     
     ''' 
     l = 0.5*(max(dx)-min(dx))
     plt.figure()
     plt.scatter(dx,dy,c=np.log10(wl))
     plt.plot(0,0,'ro')
     plt.plot([0.,l*np.cos(theta[-1])],[0,l*np.sin(theta[-1])],label = 'All members')
     plt.plot([0.,l*np.cos(theta_lum[-1])],[0,l*np.sin(theta_lum[-1])],label = 'Weighted by luminosity')
     plt.plot([0.,l*np.cos(theta_pcutdl[-1])],[0,l*np.sin(theta_pcutdl[-1])],label = 'Weighted by p_membership')
     plt.plot([0.,l*np.cos(theta_pcut[-1])],[0,l*np.sin(theta_pcut[-1])],label = 'p_membership > 0.5')
     plt.plot([0.,l*np.cos(theta_pcutl[-1])],[0,l*np.sin(theta_pcutl[-1])],label = 'p_membership > 0.5 + wL')
     plt.plot([0.,l*np.cos(theta_pcutd[-1])],[0,l*np.sin(theta_pcutd[-1])],label = 'p_membership > 0.5 + wL')
     # plt.legend()
     plt.savefig(folder+str(j)+'.png',format='png',bbox_inches='tight')
      '''

X = np.array(X)
Y = np.array(Y)

'''
tbhdu = fits.BinTableHDU.from_columns(
        [fits.Column(name='e', format='D', array=e),
        fits.Column(name='theta', format='D', array=theta),
        fits.Column(name='e_wlum', format='D', array=e_lum),
        fits.Column(name='theta_wlum', format='D', array=theta_lum),
        fits.Column(name='e_pcut_wdl', format='D', array=e_pcutdl),
        fits.Column(name='theta_pcut_wdl', format='D', array=theta_pcutdl),
        fits.Column(name='e_pcut', format='D', array=e_pcut),
        fits.Column(name='theta_pcut', format='D', array=theta_pcut),
        fits.Column(name='e_pcut_wlum', format='D', array=e_pcutl),
        fits.Column(name='theta_pcut_wlum', format='D', array=theta_pcutl),
        fits.Column(name='e_pcut_wd', format='D', array=e_pcutd),
        fits.Column(name='theta_pcut_wd', format='D', array=theta_pcutd)])
        
tbhdu.writeto(folder+'angles_redMapper.fits',overwrite=True)        


tbhdu = fits.BinTableHDU.from_columns(
        [fits.Column(name='X', format='D', array=X),
        fits.Column(name='Y', format='D', array=Y)])
        
tbhdu.writeto(folder+'redMapper_projected_member_position.fits',overwrite=True)        

'''
