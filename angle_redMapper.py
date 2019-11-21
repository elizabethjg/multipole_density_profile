import numpy as np
# from pylab import *
from astropy.io import fits
#parameters
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70., Om0=0.3, Ode0=0.7)

folder = '/home/eli/Documentos/Astronomia/posdoc/halo-elongation/redMapper/'

members  = fits.open(folder+'redmapper_dr8_public_v6.3_members.fits')
clusters = fits.open(folder+'redmapper_dr8_public_v6.3_catalog.fits')

ides  = members[1].data['ID']
ID,c = np.unique(ides,return_counts=True)

R_cen = members[1].data['R']
RA    = members[1].data['RA']
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

e_p = np.array([])
theta_p = np.array([])

e_pcut = np.array([])
theta_pcut = np.array([])

e_pcutl = np.array([])
theta_pcutl = np.array([])

e_pcutd = np.array([])
theta_pcutd = np.array([])

plots = np.linspace(1,len(ID_c),10)

# for j in plots.astype(int):

for j in range(len(ID_c)):
     print j
     mid  = ides == ID_c[j]
     N = np.append(N,mid.sum())
     ra   = RA[mid]
     dec  = DEC[mid]
     
     # all members
     mcen = R_cen[mid] == 0
     ra0  = ra[mcen]
     dec0 = dec[mcen]
     Q11  = np.sum((ra[~mcen]-ra0)**2)/np.sum(~mcen)
     Q22  = np.sum((dec[~mcen]-dec0)**2)/np.sum(~mcen)
     Q12  = np.sum((ra[~mcen]-ra0)*(dec[~mcen]-dec0))/np.sum(~mcen)
     E1 = (Q11-Q22)/(Q11+Q22)
     E2 = (2.*Q12)/(Q11+Q22)
     e = np.append(e,np.sqrt(E1**2 + E2**2))
     theta = np.append(theta,np.arctan2(E2,E1)/2.)
     
     # weighted by luminosity
     wl    = Lum_r[mid]
     mcen = R_cen[mid] == 0
     ra0  = ra[mcen]
     dec0 = dec[mcen]
     Q11  = np.sum((ra[~mcen]-ra0)**2*wl[~mcen])/np.sum(wl[~mcen])
     Q22  = np.sum((dec[~mcen]-dec0)**2*wl[~mcen])/np.sum(wl[~mcen])
     Q12  = np.sum((ra[~mcen]-ra0)*(dec[~mcen]-dec0)*wl[~mcen])/np.sum(wl[~mcen])
     E1 = (Q11-Q22)/(Q11+Q22)
     E2 = (2.*Q12)/(Q11+Q22)
     e_lum = np.append(e_lum,np.sqrt(E1**2 + E2**2))
     theta_lum = np.append(theta_lum,np.arctan2(E2,E1)/2.)
     
     # weighted by p_member
     wc   = P[mid]
     mcen = R_cen[mid] == 0
     ra0  = ra[mcen]
     dec0 = dec[mcen]
     Q11  = np.sum((ra[~mcen]-ra0)**2*wc[~mcen])/np.sum(wc[~mcen])
     Q22  = np.sum((dec[~mcen]-dec0)**2*wc[~mcen])/np.sum(wc[~mcen])
     Q12  = np.sum((ra[~mcen]-ra0)*(dec[~mcen]-dec0)*wc[~mcen])/np.sum(wc[~mcen])
     E1 = (Q11-Q22)/(Q11+Q22)
     E2 = (2.*Q12)/(Q11+Q22)
     e_p = np.append(e_p,np.sqrt(E1**2 + E2**2))
     theta_p = np.append(theta_p,np.arctan2(E2,E1)/2.)
     
     # p_cut
     pcut = P[mid] > 0.5
     mcen = R_cen[mid] == 0
     ra0  = ra[mcen]
     dec0 = dec[mcen]
     Q11  = np.sum((ra[(~mcen)*pcut]-ra0)**2)/np.sum((~mcen)*pcut)
     Q22  = np.sum((dec[(~mcen)*pcut]-dec0)**2)/np.sum((~mcen)*pcut)
     Q12  = np.sum((ra[(~mcen)*pcut]-ra0)*(dec[~mcen*pcut]-dec0))/np.sum((~mcen)*pcut)
     E1 = (Q11-Q22)/(Q11+Q22)
     E2 = (2.*Q12)/(Q11+Q22)
     e_pcut = np.append(e_pcut,np.sqrt(E1**2 + E2**2))
     theta_pcut = np.append(theta_pcut,np.arctan2(E2,E1)/2.)
     
     # p_cut weighted by lum
     pcut = P[mid] > 0.5
     mcen = R_cen[mid] == 0
     ra0  = ra[mcen]
     dec0 = dec[mcen]
     Q11  = np.sum((ra[(~mcen)*pcut]-ra0)**2*wl[(~mcen)*pcut])/np.sum(wl[(~mcen)*pcut])
     Q22  = np.sum((dec[(~mcen)*pcut]-dec0)**2*wl[(~mcen)*pcut])/np.sum(wl[(~mcen)*pcut])
     Q12  = np.sum((ra[(~mcen)*pcut]-ra0)*(dec[~mcen*pcut]-dec0)*wl[(~mcen)*pcut])/np.sum(wl[(~mcen)*pcut])
     E1 = (Q11-Q22)/(Q11+Q22)
     E2 = (2.*Q12)/(Q11+Q22)
     e_pcutl = np.append(e_pcutl,np.sqrt(E1**2 + E2**2))
     theta_pcutl = np.append(theta_pcutl,np.arctan2(E2,E1)/2.)    
     
     # p_cut weighted by distance
     wd   = 1./((ra - ra0)**2 + (dec -dec0)**2)
     pcut = P[mid] > 0.5
     mcen = R_cen[mid] == 0
     ra0  = ra[mcen]
     dec0 = dec[mcen]
     Q11  = np.sum((ra[(~mcen)*pcut]-ra0)**2*wd[(~mcen)*pcut])/np.sum(wd[(~mcen)*pcut])
     Q22  = np.sum((dec[(~mcen)*pcut]-dec0)**2*wd[(~mcen)*pcut])/np.sum(wd[(~mcen)*pcut])
     Q12  = np.sum((ra[(~mcen)*pcut]-ra0)*(dec[~mcen*pcut]-dec0)*wd[(~mcen)*pcut])/np.sum(wd[(~mcen)*pcut])
     E1 = (Q11-Q22)/(Q11+Q22)
     E2 = (2.*Q12)/(Q11+Q22)
     e_pcutd = np.append(e_pcutd,np.sqrt(E1**2 + E2**2))
     theta_pcutd = np.append(theta_pcutd,np.arctan2(E2,E1)/2.)         
     
     ''' 
     l = 0.5*(max(ra)-min(ra))
     plt.figure()
     plt.scatter(ra,dec,c=np.log10(wl))
     plt.plot(ra0,dec0,'ro')
     plt.plot([ra0,ra0+l*np.cos(theta[-1])],[dec0,dec0+l*np.sin(theta[-1])],label = 'All members')
     plt.plot([ra0,ra0+l*np.cos(theta_lum[-1])],[dec0,dec0+l*np.sin(theta_lum[-1])],label = 'Weighted by luminosity')
     plt.plot([ra0,ra0+l*np.cos(theta_p[-1])],[dec0,dec0+l*np.sin(theta_p[-1])],label = 'Weighted by p_membership')
     plt.plot([ra0,ra0+l*np.cos(theta_pcut[-1])],[dec0,dec0+l*np.sin(theta_pcut[-1])],label = 'p_membership > 0.5')
     plt.plot([ra0,ra0+l*np.cos(theta_pcutl[-1])],[dec0,dec0+l*np.sin(theta_pcutl[-1])],label = 'p_membership > 0.5 + wL')
     plt.legend()
     plt.savefig(folder+str(ID_c[j])+'.png',format='png',bbox_inches='tight')
      # '''

# '''
tbhdu = fits.BinTableHDU.from_columns(
        [fits.Column(name='e', format='D', array=e),
        fits.Column(name='theta', format='D', array=theta),
        fits.Column(name='e_wlum', format='D', array=e_lum),
        fits.Column(name='theta_wlum', format='D', array=theta_lum),
        fits.Column(name='e_wp', format='D', array=e_p),
        fits.Column(name='theta_wp', format='D', array=theta_p),
        fits.Column(name='e_pcut', format='D', array=e_pcut),
        fits.Column(name='theta_pcut', format='D', array=theta_pcut),
        fits.Column(name='e_pcut_wlum', format='D', array=e_pcutl),
        fits.Column(name='theta_pcut_wlum', format='D', array=theta_pcutl),
        fits.Column(name='e_pcut_wd', format='D', array=e_pcutd),
        fits.Column(name='theta_pcut_wd', format='D', array=theta_pcutd)])
        
tbhdu.writeto(folder+'angles_redMapper.fits',overwrite=True)        
'''
