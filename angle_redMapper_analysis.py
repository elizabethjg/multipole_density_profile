import numpy as np
# from pylab import *
from astropy.io import fits
from astropy.cosmology import LambdaCDM
import modified_corner
cosmo = LambdaCDM(H0=70., Om0=0.3, Ode0=0.7)
from matplotlib import rc
import medianas
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 12})

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

f, ax = plt.subplots(2, 4, sharex=True, sharey=True,figsize=(10,5))

# ax[0,0].xaxis.set_ticks([0.2, 0.4, 0.6, 0.8])
# ax[0,0].set_xticklabels([0.2, 0.4, 0.6, 0.8])

ax[0,0].text(0.45,500,labels[0])
ax[0,1].text(0.6,500,labels[1])
ax[0,2].text(0.6,500,labels[2])
ax[1,0].text(0.45,500,labels[3])
ax[1,1].text(0.6,500,labels[4])
ax[1,2].text(0.6,500,labels[5])
ax[1,3].text(0.6,500,labels[6])

ax[0,0].hist(inc[:,0],100,histtype='step',edgecolor='k')
ax[0,1].hist(inc[:,1],100,histtype='step',edgecolor='k')
ax[0,2].hist(inc[:,2],100,histtype='step',edgecolor='k')
ax[1,0].hist(inc[:,3],100,histtype='step',edgecolor='k')
ax[1,1].hist(inc[:,4],100,histtype='step',edgecolor='k')
ax[1,2].hist(inc[:,5],100,histtype='step',edgecolor='k')
ax[1,3].hist(inc[:,6],100,histtype='step',edgecolor='k')

ax[0,0].axvline(truths[0],color='C2')
ax[0,1].axvline(truths[1],color='C2')
ax[0,2].axvline(truths[2],color='C2')
ax[1,0].axvline(truths[3],color='C2')
ax[1,1].axvline(truths[4],color='C2')
ax[1,2].axvline(truths[5],color='C2')
ax[1,3].axvline(truths[6],color='C2')   

ax[0,0].axvline(truths[0]-std[0],c='C2',ls='--')
ax[0,1].axvline(truths[1]-std[1],c='C2',ls='--')
ax[0,2].axvline(truths[2]-std[2],c='C2',ls='--')
ax[1,0].axvline(truths[3]-std[3],c='C2',ls='--')
ax[1,1].axvline(truths[4]-std[4],c='C2',ls='--')
ax[1,2].axvline(truths[5]-std[5],c='C2',ls='--')
ax[1,3].axvline(truths[6]-std[6],c='C2',ls='--')   

ax[0,0].axvline(truths[0]+std[0],c='C2',ls='--')
ax[0,1].axvline(truths[1]+std[1],c='C2',ls='--')
ax[0,2].axvline(truths[2]+std[2],c='C2',ls='--')
ax[1,0].axvline(truths[3]+std[3],c='C2',ls='--')
ax[1,1].axvline(truths[4]+std[4],c='C2',ls='--')
ax[1,2].axvline(truths[5]+std[5],c='C2',ls='--')
ax[1,3].axvline(truths[6]+std[6],c='C2',ls='--')    
ax[0,3].axis('off')

ax[1,0].set_xlabel('$\epsilon$')
ax[1,1].set_xlabel('$\epsilon$')
ax[1,2].set_xlabel('$\epsilon$')
ax[1,3].set_xlabel('$\epsilon$')

ax[1,0].set_ylabel('$N$')
ax[0,0].set_ylabel('$N$')
plt.axis([0,0.89,0,750])
f.subplots_adjust(hspace=0,wspace=0)
plt.savefig(folder+'e_comparison.eps',format='eps',bbox_inches='tight')


style = ['C3','C4','C5','C3--','C4--','C5--','C6--']

f, ax = plt.subplots(figsize=(4.5,5))
f2, ax2 = plt.subplots(figsize=(4.5,5))
mask = Lambda < 150
for j in range(inc.shape[1]):
    x,y,d = medianas.plot_medianas_disp(Lambda[mask],inc[mask,j],plot=False)   
    ax.plot(x,y,style[j],label=labels[j])
    ax2.plot(x,d,style[j],label=labels[j])
    
ax.set_xlabel('$\lambda$')
ax.set_ylabel('$\epsilon$')
ax.axis([18,110,0.15,0.4])
ax.legend(ncol=2,fontsize=12)
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

f, ax = plt.subplots(2, 4, sharex=True, sharey=True,figsize=(10,10))
ax[0,0].hist(inc[:,0],100,histtype='step',edgecolor='k')
ax[0,1].hist(inc[:,1],100,histtype='step',edgecolor='k')
ax[0,2].hist(inc[:,2],100,histtype='step',edgecolor='k')
ax[0,3].hist(inc[:,3],100,histtype='step',edgecolor='k')
ax[1,0].hist(inc[:,4],100,histtype='step',edgecolor='k')
ax[1,1].hist(inc[:,5],100,histtype='step',edgecolor='k')
ax[1,2].hist(inc[:,6],100,histtype='step',edgecolor='k')

plt.savefig(folder+'theta_comparison.eps',format='pdf',bbox_inches='tight')

f2, ax2 = plt.subplots()
f, ax = plt.subplots()
mask = Lambda < 150
for j in range(inc.shape[1]):
    x,y,d = medianas.plot_medianas_disp(Lambda[mask],inc[mask,j],plot=False)   
    ax.plot(x,y)
    ax2.plot(x,d)
    
ax.set_xlabel('$\lambda$')
ax.set_ylabel('$\epsilon$')

# plt.savefig(folder+'e_lambda.eps',format='eps',bbox_inches='tight')


# fig = modified_corner.corner(inc,labels = labels,
                              # range = [(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max()),(inc.min(),inc.max())],
                              # plot_contours = False, truths = truths)

# plt.savefig(folder+'theta_comparison.pdf',format='pdf',bbox_inches='tight')
