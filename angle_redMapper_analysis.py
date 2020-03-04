import numpy as np
# from pylab import *
from astropy.io import fits
from astropy.cosmology import LambdaCDM
#import modified_corner
cosmo = LambdaCDM(H0=70., Om0=0.3, Ode0=0.7)
from matplotlib import rc
import medianas
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 12})

folder = '/home/eli/Documentos/Astronomia/posdoc/halo-elongation/redMapper/'
# folder = '/home/eli/Documentos/PostDoc/halo-elongation/redMapper/'

members  = fits.open(folder+'redmapper_dr8_public_v6.3_members.fits')[1].data
clusters = fits.open(folder+'redmapper_dr8_public_v6.3_catalog.fits')[1].data
angles   = fits.open(folder+'angles_redMapper_forprofile.fits')[1].data
used     = np.loadtxt(folder+'IDs_usedclusters.list')

mid         = np.in1d(clusters.ID,used)
mid_members = np.in1d(members.ID,used)
members      = members[mid_members]
clusters = clusters[mid]
angles   = angles[mid]

ides  = members['ID']
R_cen = members['R']#[mask]
ra    = members['RA']#[mask]
ra[ra > 275] = ra[ra>275] - 360.
dec   = members['DEC']#[mask]
ID,c = np.unique(ides,return_counts=True)

ID_c  = clusters['ID']
zspec = clusters['Z_SPEC']
zlambda = clusters['Z_LAMBDA']
zc = zspec
zc[zc<0] = zlambda[zc<0]
Lambda = clusters['LAMBDA']


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


# mask = (e<0.8)*(ewl<0.8)*(epwl<0.8)*(epwd<0.8)*(epwdl<0.8)
inc = np.zeros((len(e),7))

inc[:,0] = e
inc[:,1] = ewl
inc[:,2] = ewd
inc[:,3] = ep
inc[:,4] = epwl
inc[:,5] = epwd
inc[:,6] = epwdl

truths = np.average(inc,axis=0)
std    = np.std(inc,axis=0)



labels = ['$uniform$','$wL$',
          '$wd$','$uniform$',
          '$wL$' ,'$wd$', '$wdL$']

style = ['C3','C4','C5','C3--','C4--','C5--','C6--']

f, ax = plt.subplots(figsize=(4,5))
f2, ax2 = plt.subplots(figsize=(4,5))
mask = Lambda < 150
for j in range(inc.shape[1]-1):
    x,y,d = medianas.plot_medianas_disp(Lambda[mask],inc[mask,j],plot=False) 
    if j < 3:  
        ax.plot(x,y,style[j],label=labels[j])
        ax2.plot(x,d,style[j],label=labels[j])
    else:
        ax.plot(x,y,style[j])
        ax2.plot(x,d,style[j])
    
ax.set_xlabel('$\lambda$')
ax.set_ylabel('$\epsilon$')
ax.axis([18,110,0.15,0.4])
ax.legend(fontsize=12)
plt.savefig(folder+'e_lambda.eps',format='eps',bbox_inches='tight')

ax2.set_xlabel('$\lambda$')
ax2.set_ylabel('$\sigma_\epsilon$')
ax2.axis([18,110,0.08,0.19])
ax2.legend(ncol=2,fontsize=12)
plt.savefig(folder+'disp_e_lambda.eps',format='eps',bbox_inches='tight')

'''
fig = modified_corner.corner(inc,labels = labels,
                              range = [(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max())],
                              plot_contours = False, truths = truths)
'''



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

plt.hist(inc[:,1]-inc[:,0],20,color='C4',histtype='step')
plt.hist(inc[:,2]-inc[:,0],20,color='C5',histtype='step')
plt.hist(inc[:,3]-inc[:,0],20,color='C3',histtype='step',ls='--')
plt.hist(inc[:,4]-inc[:,0],20,color='C4',histtype='step',ls='--')
plt.hist(inc[:,5]-inc[:,0],20,color='C5',histtype='step',ls='--')

print np.mean(inc[:,0]-inc[:,2]),np.std(inc[:,0]-inc[:,2])
print np.mean(inc[:,1]-inc[:,2]),np.std(inc[:,1]-inc[:,2])
print np.mean(inc[:,3]-inc[:,2]),np.std(inc[:,3]-inc[:,2])
print np.mean(inc[:,4]-inc[:,2]),np.std(inc[:,4]-inc[:,2])
print np.mean(inc[:,5]-inc[:,2]),np.std(inc[:,5]-inc[:,2])

f2, ax2 = plt.subplots()
f, ax = plt.subplots()

for j in range(inc.shape[1]):
    x,y,d = medianas.plot_medianas_disp(Lambda,inc[:,j]-inc[:,0],plot=False)   
    ax.plot(x,y,style[j],label=labels[j])
    ax2.plot(x,d,style[j],label=labels[j])
    
ax.set_xlabel('$\lambda$')
ax.set_ylabel('$\epsilon$')
