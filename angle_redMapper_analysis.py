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




mask = (e<0.8)*(ewl<0.8)*(ewl<0.8)*(epwl<0.8)*(epwdl<0.8)
inc = np.zeros((len(e),5))

inc[:,0] = e
inc[:,1] = ewl
inc[:,2] = epwl
inc[:,3] = epwd
inc[:,4] = epwdl

truths = np.average(inc,axis=0)

inc = inc[mask]

labels = ['$\epsilon^{(1)}_{sat}$','$\epsilon^{(2)}_{sat}$',
          '$\epsilon^{(3)}_{sat}$','$\epsilon^{(4)}_{sat}$',
          '$\epsilon^{(5)}_{sat}$']

fig = modified_corner.corner(inc,labels = labels,
                              range = [(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max())],
                              plot_contours = False, truths = truths)

plt.savefig(folder+'e_comparison.pdf',format='pdf',bbox_inches='tight')

labels = [r'$\theta^{(1)}_{sat}$',r'$\theta^{(2)}_{sat}$',
          r'$\theta^{(3)}_{sat}$',r'$\theta^{(4)}_{sat}$',
          r'$\theta^{(5)}_{sat}$']

inc = np.zeros((len(e),5))
inc[:,0] = np.rad2deg(np.abs(t))
inc[:,1] = np.rad2deg(np.abs(twl))
inc[:,2] = np.rad2deg(np.abs(tpwl))
inc[:,3] = np.rad2deg(np.abs(tpwd))
inc[:,4] = np.rad2deg(np.abs(tpwdl))

truths = np.average(inc,axis=0)


fig = modified_corner.corner(inc,labels = labels,
                              range = [(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max())],
                              plot_contours = False, truths = truths)

plt.savefig(folder+'theta_comparison.pdf',format='pdf',bbox_inches='tight')
