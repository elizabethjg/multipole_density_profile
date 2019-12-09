import numpy as np
# from pylab import *
from astropy.io import fits
from astropy.cosmology import LambdaCDM
import modified_corner
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
twd    = angles['theta_wd']
tp    = angles['theta_pcut']
tpwl  = angles['theta_pcut_wlum']
tpwd  = angles['theta_pcut_wd']
tpwdl = angles['theta_pcut_wdl']

e     = angles['e']
ewl   = angles['e_wlum']
ewd    = angles['e_wd']
ep    = angles['e_pcut']
epwl  = angles['e_pcut_wlum']
epwd  = angles['e_pcut_wd']
epwdl = angles['e_pcut_wdl']




mask = (e<0.8)*(ewl<0.8)*(epwl<0.8)*(epwd<0.8)*(epwdl<0.8)
inc = np.zeros((len(e),7))

inc[:,0] = e
inc[:,1] = ewl
inc[:,2] = ewd
inc[:,3] = ep
inc[:,4] = epwl
inc[:,5] = epwd
inc[:,6] = epwdl

truths = np.average(inc,axis=0)

inc = inc[mask]

labels = ['$\epsilon^{(1)}_{sat}$','$\epsilon^{(2)}_{sat}$',
          '$\epsilon^{(3)}_{sat}$','$\epsilon^{(4)}_{sat}$',
          '$\epsilon^{(5)}_{sat}$','$\epsilon^{(6)}_{sat}$',
          '$\epsilon^{(7)}_{sat}$']

zfig = modified_corner.corner(inc,labels = labels,
                              range = [(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max())],
                              plot_contours = False, truths = truths)

plt.savefig(folder+'e_comparison.pdf',format='pdf',bbox_inches='tight')

labels = [r'$\theta^{(1)}_{sat}$',r'$\theta^{(2)}_{sat}$',
          r'$\theta^{(3)}_{sat}$',r'$\theta^{(4)}_{sat}$',
          r'$\theta^{(5)}_{sat}$',r'$\theta^{(6)}_{sat}$',
          r'$\theta^{(7)}_{sat}$']

inc = np.zeros((len(e),7))
inc[:,0] = np.rad2deg(np.abs(t))
inc[:,1] = np.rad2deg(np.abs(twl))
inc[:,2] = np.rad2deg(np.abs(twd))
inc[:,3] = np.rad2deg(np.abs(tp))
inc[:,4] = np.rad2deg(np.abs(tpwl))
inc[:,5] = np.rad2deg(np.abs(tpwd))
inc[:,6] = np.rad2deg(np.abs(tpwdl))

truths = np.average(inc,axis=0)


fig = modified_corner.corner(inc,labels = labels,
                              range = [(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max())],
                              plot_contours = False, truths = truths)

plt.savefig(folder+'theta_comparison.pdf',format='pdf',bbox_inches='tight')
