import numpy as np
# from pylab import *
from astropy.io import fits
sys.path.append('/home/eli/python_codes')
from astropy.cosmology import LambdaCDM
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 12})
import medianas
cosmo = LambdaCDM(H0=70., Om0=0.3, Ode0=0.7)

folder = '/home/eli/Documentos/Astronomia/posdoc/halo-elongation/redMapper/'
# folder = '/home/eli/Documentos/PostDoc/halo-elongation/redMapper/'

members       = fits.open(folder+'redmapper_dr8_public_v6.3_members.fits')[1].data
clusters      = fits.open(folder+'redmapper_dr8_public_v6.3_catalog.fits')[1].data
angles        = fits.open(folder+'angles_redMapper.fits')[1].data
coordinates   = fits.open(folder+'redMapper_projected_member_position.fits')[1].data
used          = np.loadtxt(folder+'IDs_usedclusters.list')

# --------------------------
# change angle for profile
# --------------------------

at = angles.theta
angles.theta[at<0.] = np.pi + at[at<0.]

at = angles.theta_wlum
angles.theta_wlum[at<0.] = np.pi + at[at<0.]

at = angles.theta_wd
angles.theta_wd[at<0.] = np.pi + at[at<0.]

at = angles.theta_pcut
angles.theta_pcut[at<0.] = np.pi + at[at<0.]

at = angles.theta_pcut_wlum
angles.theta_pcut_wlum[at<0.] = np.pi + at[at<0.]

at = angles.theta_pcut_wd
angles.theta_pcut_wd[at<0.] = np.pi + at[at<0.]

at = angles.theta_pcut_wdl
angles.theta_pcut_wdl[at<0.] = np.pi + at[at<0.]

hdu = fits.BinTableHDU(angles)
hdu.writeto(folder+'angles_redMapper_forprofile.fits',overwrite=True)


# --------------------------
angles        = fits.open(folder+'angles_redMapper.fits')[1].data
mid         = np.in1d(clusters.ID,used)
mid_members = np.in1d(members.ID,used)

members      = members[mid_members]
coordinates  = coordinates[mid_members]
clusters = clusters[mid]
angles   = angles[mid]

ides  = members.ID
R_cen = members.R
RA    = members.RA
RA[RA > 275] = RA[RA>275] - 360.
DEC   = members.DEC
P     = members.P

ID,c = np.unique(ides,return_counts=True)

ID_c  = clusters.ID
zspec = clusters.Z_SPEC
zlambda = clusters.Z_LAMBDA
zc = zspec
zc[zc<0] = zlambda[zc<0]
Lambda = np.repeat(clusters.LAMBDA,c)
R_lambda = ((Lambda/100.)**(0.2))/0.7


#######################################
# HISTOGRAMS redMaPPer sample
f, ax = plt.subplots(2, 1, figsize=(4,6))
ax[0].hist(zc,50,histtype='step',color='C5')
ax[0].set_xlabel(r'$z_{cluster}$',fontsize = 14)
ax[0].set_ylabel(r'$N$')
ax[0].axis([0.1,0.4,0.,100])
ax[0].axvline(0.313,c='C4')
ax[1].hist(clusters.LAMBDA,50,histtype='step',color='C5')
ax[1].set_xlabel(r'$\lambda$',fontsize = 14)
ax[1].axvline(27.982,c='C4')
ax[1].set_ylabel(r'$N$')
plt.subplots_adjust(hspace = 0.25)
ax[1].axis([20,150,0.,400])
# f.subplots_adjust(hspace=0,wspace=0)
plt.savefig(folder+'hist.eps',format='eps',bbox_inches='tight')
#######################################

D_ang    = np.array(cosmo.angular_diameter_distance(zc))
kpcscale = D_ang*(((1.0/3600.0)*np.pi)/180.0)*1000.0
KPCSCALE = np.repeat(kpcscale,c)

zc = np.repeat(zc,c)
D_lum    = np.array(cosmo.luminosity_distance(zc))*1.e6
MAG_abs  = members.MODEL_MAG_R + 5.0-5.0*np.log10(D_lum)
Lum_r = 10.**(-0.4*MAG_abs)

c_gr = members.MODEL_MAG_G - members.MODEL_MAG_R
c_ri = members.MODEL_MAG_R - members.MODEL_MAG_I

mcen = R_cen == 0.



RA0  = np.repeat(RA[mcen],c)
DEC0 = np.repeat(DEC[mcen],c)

t     = np.repeat(angles['theta'],c)
twl   = np.repeat(angles['theta_wlum'],c)
twd   = np.repeat(angles['theta_wd'],c)
tp    = np.repeat(angles['theta_pcut'],c)
tpwl  = np.repeat(angles['theta_pcut_wlum'],c)
tpwd  = np.repeat(angles['theta_pcut_wd'],c)
tpwdl = np.repeat(angles['theta_pcut_wdl'],c)


    
e     = np.repeat(angles['e'],c)
ewl   = np.repeat(angles['e_wlum'],c)
ewd   = np.repeat(angles['e_wd'],c)
ep    = np.repeat(angles['e_pcut'],c)
epwl  = np.repeat(angles['e_pcut_wlum'],c)
epwd  = np.repeat(angles['e_pcut_wd'],c)
epwdl = np.repeat(angles['e_pcut_wdl'],c)
# '''

dx = (coordinates['X']*KPCSCALE*1.e-3)/R_lambda
dy = (coordinates['Y']*KPCSCALE*1.e-3)/R_lambda

mask = np.isfinite(dx)
dx = dx[mask]
dy = dy[mask]

x_t  = (dx*np.cos(t[mask]) + dy*np.sin(t[mask])) 
y_t = (-1.*dx*np.sin(t[mask]) + dy*np.cos(t[mask])) 

x_l  = (dx*np.cos(twl[mask]) + dy*np.sin(twl[mask])) 
y_l = (-1.*dx*np.sin(twl[mask]) + dy*np.cos(twl[mask])) 

x_d  = (dx*np.cos(twd[mask]) + dy*np.sin(twd[mask])) 
y_d = (-1.*dx*np.sin(twd[mask]) + dy*np.cos(twd[mask])) 

x_p  = (dx*np.cos(tp[mask]) + dy*np.sin(tp[mask])) 
y_p = (-1.*dx*np.sin(tp[mask]) + dy*np.cos(tp[mask])) 

x_pwl  = (dx*np.cos(tpwl[mask]) + dy*np.sin(tpwl[mask])) 
y_pwl = (-1.*dx*np.sin(tpwl[mask]) + dy*np.cos(tpwl[mask])) 

x_pwd  = (dx*np.cos(tpwd[mask]) + dy*np.sin(tpwd[mask])) 
y_pwd = (-1.*dx*np.sin(tpwd[mask]) + dy*np.cos(tpwd[mask])) 

x_pwdl  = (dx*np.cos(tpwdl[mask]) + dy*np.sin(tpwdl[mask])) 
y_pwdl = (-1.*dx*np.sin(tpwdl[mask]) + dy*np.cos(tpwdl[mask])) 


lgrid = 20
Max = 1800


'''
f, ax = plt.subplots(2, 4, sharex=True, sharey=True,figsize=(10,10))
ax[0,0].hexbin(dx,dy,gridsize=lgrid,extent=(-0.4,0.4,-0.4,0.4),vmax=Max,cmap='pink')
ax[0,1].hexbin(x_t,y_t,gridsize=lgrid,extent=(-0.4,0.4,-0.4,0.4),vmax=Max,cmap='pink')
ax[0,2].hexbin(x_l,y_l,gridsize=lgrid,extent=(-0.4,0.4,-0.4,0.4),vmax=Max,cmap='pink')
ax[0,3].hexbin(x_d,y_d,gridsize=lgrid,extent=(-0.4,0.4,-0.4,0.4),vmax=Max,cmap='pink')
ax[1,0].hexbin(x_p,y_p,gridsize=lgrid,extent=(-0.4,0.4,-0.4,0.4),vmax=Max,cmap='pink')
ax[1,1].hexbin(x_pwl,y_pwl,gridsize=lgrid,extent=(-0.4,0.4,-0.4,0.4),vmax=Max,cmap='pink')
ax[1,2].hexbin(x_pwd,y_pwd,gridsize=lgrid,extent=(-0.4,0.4,-0.4,0.4),vmax=Max,cmap='pink')
ax[1,3].hexbin(x_pwdl,y_pwdl,gridsize=lgrid,extent=(-0.4,0.4,-0.4,0.4),vmax=Max,cmap='pink')
f.subplots_adjust(hspace=0,wspace=0)
'''

xedges = np.linspace(-2.,2.,50)
yedges = np.linspace(-2.,2.,50)
xcenters = (xedges[:-1] + xedges[1:]) / 2.
ycenters = (yedges[:-1] + yedges[1:]) / 2.
X,Y = np.meshgrid(xcenters,ycenters)
f, ax = plt.subplots(2, 4, sharey=False,figsize=(8,4))


levels = np.linspace(150,1000,10)
H, xedges, yedges = np.histogram2d(dx, dy, bins=(xedges, yedges))#,weights=np.log10(Lum_r))
ax[0,0].contour(X, Y, H.T,levels,cmap='plasma')
ax[0,0].set_ylabel('$x_2/R_{\lambda}$',fontsize = '14')
ax[0,0].set_xlabel('$x_1/R_{\lambda}$',fontsize = '14')
ax[0,0].axis([-1.2,1.2,-1.2,1.2])

H, xedges, yedges = np.histogram2d(x_t,y_t, bins=(xedges, yedges))#,weights=np.log10(Lum_r))
ax[0,1].contour(X, Y, H.T,levels,cmap='plasma')
ax[0,1].axis([-1.2,1.2,-1.2,1.2])
ax[0,1].text(0.7,0.8,'$\phi_1$',fontsize = 14)
plt.setp(ax[0,1].get_xticklabels(), visible=False)
plt.setp(ax[0,1].get_yticklabels(), visible=False)

H, xedges, yedges = np.histogram2d(x_l,y_l, bins=(xedges, yedges))#,weights=np.log10(Lum_r))
ax[0,2].contour(X, Y, H.T,levels,cmap='plasma')
ax[0,2].axis([-1.2,1.2,-1.2,1.2])
ax[0,2].text(0.7,0.8,'$\phi_L$',fontsize = 14)
plt.setp(ax[0,2].get_xticklabels(), visible=False)
plt.setp(ax[0,2].get_yticklabels(), visible=False)

H, xedges, yedges = np.histogram2d(x_d,y_d, bins=(xedges, yedges))#,weights=np.log10(Lum_r))
ax[0,3].contour(X, Y, H.T,levels,cmap='plasma')
ax[0,3].axis([-1.2,1.2,-1.2,1.2])
ax[0,3].text(0.7,0.8,r'$\phi_d$',fontsize = 14)
plt.setp(ax[0,3].get_xticklabels(), visible=False)
plt.setp(ax[0,3].get_yticklabels(), visible=False)

ax[1,0].axis('off')

H, xedges, yedges = np.histogram2d(x_p,y_p, bins=(xedges, yedges))#,weights=np.log10(Lum_r))
ax[1,1].contour(X, Y, H.T,levels,cmap='plasma')
ax[1,1].axis([-1.2,1.2,-1.2,1.2])
ax[1,1].text(0.7,0.8,u'$\phi^*_1$',fontsize = 14)
ax[1,1].set_ylabel('$x_2/R_{\lambda}$',fontsize = '14')
ax[1,1].set_xlabel('$x_1/R_{\lambda}$',fontsize = '14')
ax[1,1].yaxis.set_ticks([-1,0])

H, xedges, yedges = np.histogram2d(x_pwl,y_pwl, bins=(xedges, yedges))#,weights=np.log10(Lum_r))
ax[1,2].contour(X, Y, H.T,levels,cmap='plasma')
ax[1,2].axis([-1.2,1.2,-1.2,1.2])
ax[1,2].text(0.7,0.8,u'$\phi^*_L$',fontsize = 14)
ax[1,2].set_xlabel('$x_1/R_{\lambda}$',fontsize = '14')
plt.setp(ax[1,2].get_yticklabels(), visible=False)

H, xedges, yedges = np.histogram2d(x_pwd,y_pwd, bins=(xedges, yedges))#,weights=np.log10(Lum_r))
ax[1,3].contour(X, Y, H.T,levels,cmap='plasma')
ax[1,3].axis([-1.2,1.2,-1.2,1.2])
ax[1,3].text(0.7,0.8,u'$\phi^*_d$',fontsize = 14)
ax[1,3].set_xlabel('$x_1/R_{\lambda}$',fontsize = '14')
plt.setp(ax[1,3].get_yticklabels(), visible=False)

f.subplots_adjust(hspace=0,wspace=0)
plt.savefig(folder+'contours.pdf',format='pdf',bbox_inches='tight')

# '''
#########################
#  DISTANCE DISTRIBUTION

Rscaled = (R_cen/0.7)/R_lambda
f, ax = plt.subplots(1, 1, figsize=(5,4),sharex=True)

ax.hist(Rscaled[~mcen],np.linspace(0,1.2,50),histtype='step',density=True,color='C3',)
ax.axvline(np.average(Rscaled[~mcen]),color='C3',label=u'$\sqrt{w_1}$')
ax.hist(Rscaled[~mcen],np.linspace(0,1.2,50),histtype='step',weights=np.sqrt(Lum_r[~mcen]),density=True,color='C4',) 
ax.axvline(np.average(Rscaled[~mcen],weights=np.sqrt(Lum_r[~mcen])),color='C4',label=u'$\sqrt{w_L}$')
ax.hist(Rscaled[~mcen],np.linspace(0,1.2,50),histtype='step',weights=np.sqrt(1./(dx**2 + dy**2))[~mcen],density=True,color='C5',) 
ax.axvline(np.average(Rscaled[~mcen],weights=np.sqrt(1./(dx**2 + dy**2))[~mcen]),color='C5',label=u'$\sqrt{w_d}$')
ax.hist(Rscaled[(P>0.5)*(~mcen)],np.linspace(0,1.2,50),histtype='step',density=True,linestyle='dashed',color='C3',) 
ax.axvline(np.average(Rscaled[(P>0.5)*(~mcen)]),color='C3',ls='--')
ax.hist(Rscaled[(P>0.5)*(~mcen)],np.linspace(0,1.2,50),histtype='step',weights=np.sqrt(Lum_r[(P>0.5)*(~mcen)]),density=True,linestyle='dashed',color='C4',) 
ax.axvline(np.average(Rscaled[(P>0.5)*(~mcen)],weights=np.sqrt(Lum_r[(P>0.5)*(~mcen)])),color='C4',ls='--')
ax.hist(Rscaled[(P>0.5)*(~mcen)],np.linspace(0,1.2,50),histtype='step',weights=np.sqrt(1./(dx**2 + dy**2))[(P>0.5)*(~mcen)],density=True,linestyle='dashed',color='C5',) 
ax.axvline(np.average(Rscaled[(P>0.5)*(~mcen)],weights=np.sqrt(1./(dx**2 + dy**2))[(P>0.5)*(~mcen)]),color='C5',ls='--')
ax.set_ylabel('$n$',fontsize = '14')
ax.legend(frameon=False)
ax.axis([0.,1.2,0,4.7])
ax.set_xlabel('$r/R_{\lambda}$',fontsize = '14')
ax.set_ylabel('$n$',fontsize = '14')

plt.savefig(folder+'dist.eps',format='eps',bbox_inches='tight')

#########################
# '''

'''
levels = np.linspace(H.min(),10000.,15)

f, ax = plt.subplots(2, 4, sharex=True, sharey=True,figsize=(8,4))
H, xedges, yedges = np.histogram2d(dx, dy, bins=(xedges, yedges))
ax[0,0].contour(X, Y, H.T,levels)
ax[0,0].set_ylabel('$x_2/R_{\lambda}$',fontsize = '14')

H, xedges, yedges = np.histogram2d(x_t,y_t, bins=(xedges, yedges))
ax[0,1].contour(X, Y, H.T,levels)

H, xedges, yedges = np.histogram2d(x_l,y_l, bins=(xedges, yedges))
ax[0,2].contour(X, Y, H.T,levels)

H, xedges, yedges = np.histogram2d(x_d,y_d, bins=(xedges, yedges))
ax[0,3].contour(X, Y, H.T,levels)

H, xedges, yedges = np.histogram2d(x_p,y_p, bins=(xedges, yedges))
ax[1,0].contour(X, Y, H.T,levels)
ax[1,0].set_xlabel('$x_1/R_{\lambda}$',fontsize = '14')
ax[1,0].set_ylabel('$x_2/R_{\lambda}$',fontsize = '14')

H, xedges, yedges = np.histogram2d(x_pwl,y_pwl, bins=(xedges, yedges))
ax[1,1].contour(X, Y, H.T,levels)
ax[1,1].set_xlabel('$x_1/R_{\lambda}$',fontsize = '14')

H, xedges, yedges = np.histogram2d(x_pwd,y_pwd, bins=(xedges, yedges))
ax[1,2].contour(X, Y, H.T,levels)
ax[1,2].set_xlabel('$x_1/R_{\lambda}$',fontsize = '14')

H, xedges, yedges = np.histogram2d(x_pwdl,y_pwdl, bins=(xedges, yedges))
ax[1,3].contour(X, Y, H.T,levels)
ax[1,3].set_xlabel('$x_1/R_{\lambda}$',fontsize = '14')

f.subplots_adjust(hspace=0,wspace=0)
plt.savefig(folder+'contours.eps',format='eps',bbox_inches='tight')

'''


