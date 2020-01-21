import numpy as np
# from pylab import *
from astropy.io import fits
#parameters
from astropy.cosmology import LambdaCDM
from maria_func import *
from make_profile import *
cosmo = LambdaCDM(H0=70., Om0=0.3, Ode0=0.7)
#parameters
cvel = 299792458;   # Speed of light (m.s-1)
G    = 6.670e-11;   # Gravitational constant (m3.kg-1.s-2)
pc   = 3.085678e16; # 1 pc (m)
Msun = 1.989e30 # Solar mass (kg)
'''
z_back   = 0.1
odds_min = 0.5
zmin     = 0.1
zmax     = 0.33
RIN      = 100.
ROUT     = 10000.
ndots    = 20.
sample   = 'bin4'
lmin     = 39.7
lmax     = 145.
'''

def profile_redMapper(sample,lmin,lmax,zmin = 0.1, zmax = 0.33,
                      z_back = 0.1, odds_min = 0.5,
                      RIN = 100., ROUT = 10000., ndots = 20.):

     folder = '/mnt/clemente/lensing/redMaPPer/'
     
     cfht     = fits.open(folder+'gx_CFHT_redMapper.fits')[1].data
     kids     = fits.open(folder+'gx_KiDS_redMapper.fits')[1].data
     cs82     = fits.open(folder+'gx_CS82_redMapper.fits')[1].data
     angles   = fits.open(folder+'angles_redMapper.fits')[1].data
     clusters = fits.open(folder+'redmapper_dr8_public_v6.3_catalog.fits')[1].data
     borderid = np.loadtxt(folder+'redMapperID_border.list')
     
     mkids = ~np.in1d(kids.ID,cfht.ID)
     mcs82 = ~np.in1d(cs82.ID,cfht.ID)
     
     kids = kids[mkids]
     cs82 = cs82[mcs82]
     
     # MATCH ANGLES
     
     ID = np.concatenate((cfht.ID,kids.ID,cs82.ID))
     IDc = clusters.ID
     
     ides,index,c = np.unique(ID,return_index=True,return_counts=True)
     angles_in = angles[np.in1d(IDc,ID)]
     sindex = np.argsort(index)
     angles = np.repeat(angles_in[sindex],c[sindex])
     
     #-----------------------------------------------
     
     lamb = np.concatenate((cfht.LAMBDA,kids.LAMBDA,cs82.LAMBDA))
     
     Z_B  = np.concatenate((cfht.Z_B,kids.Z_B,cs82.Z_B))
     ODDS = np.concatenate((cfht.ODDS,kids.ODDS,cs82.ODDS))
     zlambda = np.concatenate((cfht.Z_LAMBDA,kids.Z_LAMBDA,cs82.Z_LAMBDA))
     zspec   = np.concatenate((cfht.Z_SPEC,kids.Z_SPEC,cs82.Z_SPEC))
     Z_c      = zspec
     Z_c[Z_c<0] = zlambda[Z_c<0]
     
     
     mask_back = (Z_B > (Z_c + z_back))*(ODDS >= odds_min)#*(Z_B > (Z_c + s95/2.))
     mask_lens = (lamb >= lmin)*(lamb < lmax)*(Z_c >= zmin)*(Z_c < zmax)*(~np.in1d(ID,borderid))
     mask = mask_back*mask_lens
     
     Nclusters = len(np.unique(ID[mask]))
     
     ra     = np.concatenate((cfht.RAJ2000,kids.RAJ2000,cs82.RAJ2000))[mask]
     dec    = np.concatenate((cfht.DECJ2000,kids.DECJ2000,cs82.DECJ2000))[mask]
     Z_c    = Z_c[mask]
          
     e1     = np.concatenate((cfht.e1,kids.e1,cs82.e1))[mask]
     e2     = np.concatenate((cfht.e2,kids.e2,cs82.e2))[mask]
     
          
     ALFA0  = np.concatenate((cfht.RA,kids.RA,cs82.RA))[mask]
     DELTA0 = np.concatenate((cfht.DEC,kids.DEC,cs82.DEC))[mask]
     
     dls  = np.concatenate((cfht.DLS,kids.DLS,cs82.DLS))[mask]
     ds   = np.concatenate((cfht.DS,kids.DS,cs82.DS))[mask]
     dl   = np.concatenate((cfht.DL,kids.DL,cs82.DL))[mask]
     
     peso = np.concatenate((cfht.weight,kids.weight,cs82.weight))[mask]
     
     m    = np.concatenate((cfht.m,kids.m,cs82.m))[mask]
     
     KPCSCALE   = dl*(((1.0/3600.0)*np.pi)/180.0)*1000.0
     BETA_array = dls/ds
     beta       = BETA_array.mean()
     
     Dl = dl*1.e6*pc
     sigma_c = (((cvel**2.0)/(4.0*np.pi*G*Dl))*(1./BETA_array))*(pc**2/Msun)
     peso=peso/(sigma_c**2)
     
     print 'BETA.mean',beta
     
     SIGMAC = (((cvel**2.0)/(4.0*np.pi*G*Dl.mean())))*(pc**2/Msun)
     
     print 'SIGMA_C', SIGMAC
     
     
     rads, theta, test1,test2 = eq2p2(np.deg2rad(ra),
                              np.deg2rad(dec),
                              np.deg2rad(ALFA0),
                              np.deg2rad(DELTA0))
     
     
     #Correct polar angle for e1, e2
     theta = theta+np.pi/2.
     
     #get tangential ellipticities 
     et = (-e1*np.cos(2*theta)-e2*np.sin(2*theta))*sigma_c
     #get cross ellipticities
     ex = (-e1*np.sin(2*theta)+e2*np.cos(2*theta))*sigma_c
     
     
     r=np.rad2deg(rads)*3600*KPCSCALE
     
     
     zmean    = (Z_c).mean()
     zdisp    = (Z_c).std()
     H        = cosmo.H(zmean).value/(1.0e3*pc) #H at z_pair s-1 
     roc      = (3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_pair (kg.m-3)
     roc_mpc  = roc*((pc*1.0e6)**3.0)
     D_ang    = cosmo.angular_diameter_distance(zmean)
     kpcscale = D_ang*(((1.0/3600.0)*np.pi)/180.0)*1000.0
     
     
     print '---------------------------------------------------------'
     print '             COMPUTING THE SHEAR PROFILES                '
     print '========================================================='
     
     
     profile = shear_profile_log(RIN,ROUT,r,et,ex,peso,m,sigma_c,ndots,
                              booterror_flag=True)
                              
     print 'Now is trying to fit a NFW profile...'
     
     try:
          nfw          = NFW_stack_fit(profile[0]/1.e3,profile[1],profile[5],zmean,roc)
     except:
          nfw          = [-999.,0.,-100.,[0.,0.],[0.,0.],-999.,0.]
                              
     
     M200_NFW   = (800.0*np.pi*roc_mpc*(nfw[0]**3))/(3.0*Msun)
     
     f1=open('profile_'+sample+'.cat','w')
     f1.write('# Nclusters = '+str('%8i' % Nclusters)+' \n')
     f1.write('# M200 = '+str('%.2f' % (M200_NFW/1.e14))+' \n')
     f1.write('# z_mean = '+str('%.2f' % zmean)+' \n')
     f1.write('# z_back = '+str('%.2f' % z_back)+' \n')
     f1.write('# odds_min = '+str('%.1f' % odds_min)+' \n')
     f1.write('# l_min = '+str('%.1f' % lmin)+' \n')
     f1.write('# l_max = '+str('%.1f' % lmax)+' \n')
     f1.write('# z_min = '+str('%.1f' % zmin)+' \n')
     f1.write('# z_max = '+str('%.1f' % zmax)+' \n')
     f1.write('# R,shear,err_et,cero,err_ex \n')
     profile = np.column_stack((profile[0]*1.0e-3,profile[1],profile[5],profile[2],profile[6]))
     np.savetxt(f1,profile,fmt = ['%12.6f']*5)
     f1.close()
     
     def write_profile(angle,sample):
          
          print sample
     
          profile = quadrupole_profile_log(RIN,ROUT,r,et,ex,
                                   peso,m,sigma_c,angle,
                                   ndots,booterror_flag=True)
          
          f1=open('profile_'+sample+'.cat','w')
          f1.write('# Nclusters = '+str('%8i' % Nclusters)+' \n')
          f1.write('# M200 = '+str('%.2f' % (M200_NFW/1.e14))+' \n')
          f1.write('# z_mean = '+str('%.2f' % zmean)+' \n')
          f1.write('# z_back = '+str('%.2f' % z_back)+' \n')
          f1.write('# odds_min = '+str('%.1f' % odds_min)+' \n')
          f1.write('# l_min = '+str('%.1f' % lmin)+' \n')
          f1.write('# l_max = '+str('%.1f' % lmax)+' \n')
          f1.write('# z_min = '+str('%.1f' % zmin)+' \n')
          f1.write('# z_max = '+str('%.1f' % zmax)+' \n')
          f1.write('# R,shear,err_et,cero,err_ex \n')
          profile = np.column_stack((profile[0]*1.0e-3,profile[1],profile[5],profile[2],profile[6]))
          np.savetxt(f1,profile,fmt = ['%12.6f']*5)
          f1.close()
     
     write_profile(angles.theta[mask]+np.pi/2.,sample+'_t')
     write_profile(angles.theta_wlum[mask]+np.pi/2.,sample+'_twl')
     write_profile(angles.theta_wd[mask]+np.pi/2.,sample+'_twd')
     write_profile(angles.theta_pcut[mask]+np.pi/2.,sample+'_tp')
     write_profile(angles.theta_pcut_wlum[mask]+np.pi/2.,sample+'_tpwl')
     write_profile(angles.theta_pcut_wd[mask]+np.pi/2.,sample+'_tpwd')
     write_profile(angles.theta_pcut_wdl[mask]+np.pi/2.,sample+'_tpwdl')

def profile_redMapper_indcat(name_cat,sample,lmin,lmax,zmin = 0.1, zmax = 0.33,
                      z_back = 0.1, odds_min = 0.5,
                      RIN = 100., ROUT = 10000., ndots = 20.):

     folder = '/mnt/clemente/lensing/redMaPPer/'
     
     backgx   = fits.open(folder+name_cat)[1].data
     angles   = fits.open(folder+'angles_redMapper.fits')[1].data
     clusters = fits.open(folder+'redmapper_dr8_public_v6.3_catalog.fits')[1].data
     borderid = np.loadtxt(folder+'redMapperID_border.list')

     # MATCH ANGLES
     
     ID = backgx.ID
     IDc = clusters.ID
     
     ides,index,c = np.unique(ID,return_index=True,return_counts=True)
     angles_in = angles[np.in1d(IDc,ID)]
     sindex = np.argsort(index)
     angles = np.repeat(angles_in[sindex],c[sindex])
     
     
     #-----------------------------------------------
     
     lamb = backgx.LAMBDA
     
     Z_B  = backgx.Z_B
     ODDS = backgx.ODDS
     zlambda = backgx.Z_LAMBDA
     zspec   = backgx.Z_SPEC
     Z_c      = zspec
     Z_c[Z_c<0] = zlambda[Z_c<0]
     
     
     mask_back = (Z_B > (Z_c + z_back))*(ODDS >= odds_min)#*(Z_B > (Z_c + s95/2.))
     mask_lens = (lamb >= lmin)*(lamb < lmax)*(Z_c >= zmin)*(Z_c < zmax)*(~np.in1d(ID,borderid))
     mask = mask_back*mask_lens
     
     Nclusters = len(np.unique(ID[mask]))
     
     ra     = backgx.RAJ2000[mask]
     dec    = backgx.DECJ2000[mask]
     Z_c    = Z_c[mask]
     lamb   = lamb[mask]     
     e1     = backgx.e1[mask]
     e2     = backgx.e2[mask]
     
          
     ALFA0  = backgx.RA[mask]
     DELTA0 = backgx.DEC[mask]
     
     dls  = backgx.DLS[mask]
     ds   = backgx.DS[mask]
     dl   = backgx.DL[mask]
     
     peso = backgx.weight[mask]
     m    = backgx.m[mask]
     
     KPCSCALE   = dl*(((1.0/3600.0)*np.pi)/180.0)*1000.0
     BETA_array = dls/ds
     beta       = BETA_array.mean()
     
     Dl = dl*1.e6*pc
     sigma_c = (((cvel**2.0)/(4.0*np.pi*G*Dl))*(1./BETA_array))*(pc**2/Msun)
     
     
     print 'BETA.mean',beta
     
     SIGMAC = (((cvel**2.0)/(4.0*np.pi*G*Dl.mean())))*(pc**2/Msun)
     
     print 'SIGMA_C', SIGMAC
     
     
     rads, theta, test1,test2 = eq2p2(np.deg2rad(ra),
                              np.deg2rad(dec),
                              np.deg2rad(ALFA0),
                              np.deg2rad(DELTA0))
     
     
     #Correct polar angle for e1, e2
     theta = theta+np.pi/2.
     
     #get tangential ellipticities 
     et = (-e1*np.cos(2*theta)-e2*np.sin(2*theta))*sigma_c
     #get cross ellipticities
     ex = (-e1*np.sin(2*theta)+e2*np.cos(2*theta))*sigma_c
     
     
     r=np.rad2deg(rads)*3600*KPCSCALE
     
     
     zmean    = (Z_c).mean()
     zdisp    = (Z_c).std()
     H        = cosmo.H(zmean).value/(1.0e3*pc) #H at z_pair s-1 
     roc      = (3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_pair (kg.m-3)
     roc_mpc  = roc*((pc*1.0e6)**3.0)
     D_ang    = cosmo.angular_diameter_distance(zmean)
     kpcscale = D_ang*(((1.0/3600.0)*np.pi)/180.0)*1000.0
     
     
     print '---------------------------------------------------------'
     print '             COMPUTING THE SHEAR PROFILES                '
     print '========================================================='
     
     
     profile = shear_profile_log(RIN,ROUT,r,et,ex,peso,m,sigma_c,ndots,
                              booterror_flag=True)
                              
     print 'Now is trying to fit a NFW profile...'
     
     try:
          nfw          = NFW_stack_fit(profile[0]/1.e3,profile[1],profile[5],zmean,roc)
     except:
          nfw          = [-999.,0.,-100.,[0.,0.],[0.,0.],-999.,0.]
                              
     
     M200_NFW   = (800.0*np.pi*roc_mpc*(nfw[0]**3))/(3.0*Msun)
     
     f1=open('profile_'+sample+'.cat','w')
     f1.write('# Nclusters = '+str('%8i' % Nclusters)+' \n')
     f1.write('# M200 = '+str('%.2f' % (M200_NFW/1.e14))+' \n')
     f1.write('# z_mean = '+str('%.2f' % zmean)+' \n')
     f1.write('# z_back = '+str('%.2f' % z_back)+' \n')
     f1.write('# odds_min = '+str('%.1f' % odds_min)+' \n')
     f1.write('# l_min  = '+str('%.1f' % lmin)+' \n')
     f1.write('# l_max  = '+str('%.1f' % lmax)+' \n')
     f1.write('# l_mean = '+str('%.1f' % np.mean(lamb))+' \n')
     f1.write('# z_min = '+str('%.1f' % zmin)+' \n')
     f1.write('# z_max = '+str('%.1f' % zmax)+' \n')
     f1.write('# R,shear,err_et,cero,err_ex \n')
     profile = np.column_stack((profile[0]*1.0e-3,profile[1],profile[5],profile[2],profile[6]))
     np.savetxt(f1,profile,fmt = ['%12.6f']*5)
     f1.close()
     
     def write_profile(angle,sample):
          
          print sample
     
          profile = quadrupole_profile_log(RIN,ROUT,r,et,ex,
                                   peso,m,sigma_c,angle,
                                   ndots,booterror_flag=True)
          
          f1=open('profile_'+sample+'.cat','w')
          f1.write('# Nclusters = '+str('%8i' % Nclusters)+' \n')
          f1.write('# M200 = '+str('%.2f' % (M200_NFW/1.e14))+' \n')
          f1.write('# z_mean = '+str('%.2f' % zmean)+' \n')
          f1.write('# z_back = '+str('%.2f' % z_back)+' \n')
          f1.write('# odds_min = '+str('%.1f' % odds_min)+' \n')
          f1.write('# l_min  = '+str('%.1f' % lmin)+' \n')
          f1.write('# l_max  = '+str('%.1f' % lmax)+' \n')
          f1.write('# l_mean = '+str('%.1f' % np.mean(lamb))+' \n')
          f1.write('# z_min = '+str('%.1f' % zmin)+' \n')
          f1.write('# z_max = '+str('%.1f' % zmax)+' \n')
          f1.write('# R,shear,err_et,cero,err_ex \n')
          profile = np.column_stack((profile[0]*1.0e-3,profile[1],profile[5],profile[2],profile[6]))
          np.savetxt(f1,profile,fmt = ['%12.6f']*5)
          f1.close()
     
     write_profile(angles.theta[mask]+np.pi/2.,sample+'_t')
     write_profile(angles.theta_wlum[mask]+np.pi/2.,sample+'_twl')
     write_profile(angles.theta_wd[mask]+np.pi/2.,sample+'_twd')
     write_profile(angles.theta_pcut[mask]+np.pi/2.,sample+'_tp')
     write_profile(angles.theta_pcut_wlum[mask]+np.pi/2.,sample+'_tpwl')
     write_profile(angles.theta_pcut_wd[mask]+np.pi/2.,sample+'_tpwd')
     write_profile(angles.theta_pcut_wdl[mask]+np.pi/2.,sample+'_tpwdl')

profile_redMapper('original_bin4',39.7,145.,RIN=150.,ROUT=10000.,ndots=20)
profile_redMapper('original_total',0.,145.,RIN=150.,ROUT=10000.,ndots=20)
profile_redMapper('original_bin1',20.,23.42,RIN=150.,ROUT=10000.,ndots=20)
profile_redMapper('original_bin2',23.42,28.3,RIN=150.,ROUT=10000.,ndots=20)
profile_redMapper('original_bin3',28.3,39.7,RIN=150.,ROUT=10000.,ndots=20)

profile_redMapper('bin1',22.,40.,RIN=100.,ROUT=5000.,ndots=10)
profile_redMapper('total',0.,150.,RIN=100.,ROUT=5000.,ndots=10)
profile_redMapper('bin2',40.,70.,RIN=100.,ROUT=5000.,ndots=10)
profile_redMapper('bin3',70.,150.,RIN=100.,ROUT=5000.,ndots=10)



#profile_redMapper_indcat('gx_CFHT_redMapper.fits','cfht_bin4',39.7,145.,RIN=100.,ROUT=5000.,ndots=15)
#profile_redMapper_indcat('gx_CFHT_redMapper.fits','cfht_total',0.,145.,RIN=100.,ROUT=5000.,ndots=15)
#profile_redMapper_indcat('gx_CFHT_redMapper.fits','cfht_bin1',20.,23.42,RIN=100.,ROUT=5000.,ndots=15)
#profile_redMapper_indcat('gx_CFHT_redMapper.fits','cfht_bin2',23.42,28.3,RIN=100.,ROUT=5000.,ndots=15)
#profile_redMapper_indcat('gx_CFHT_redMapper.fits','cfht_bin3',28.3,39.7,RIN=100.,ROUT=5000.,ndots=15)
                                           
#profile_redMapper_indcat('gx_KiDS_redMapper.fits','KiDS_bin4',39.7,145.,RIN=100.,ROUT=5000.,ndots=15)
#profile_redMapper_indcat('gx_KiDS_redMapper.fits','KiDS_total',0.,145.,RIN=100.,ROUT=5000.,ndots=15)
#profile_redMapper_indcat('gx_KiDS_redMapper.fits','KiDS_bin1',20.,23.42,RIN=100.,ROUT=5000.,ndots=15)
#profile_redMapper_indcat('gx_KiDS_redMapper.fits','KiDS_bin2',23.42,28.3,RIN=100.,ROUT=5000.,ndots=15)
#profile_redMapper_indcat('gx_KiDS_redMapper.fits','KiDS_bin3',28.3,39.7,RIN=100.,ROUT=5000.,ndots=15)
                                           
#profile_redMapper_indcat('gx_CS82_redMapper.fits','CS82_bin4',39.7,145.,RIN=100.,ROUT=5000.,ndots=15)
#profile_redMapper_indcat('gx_CS82_redMapper.fits','CS82_total',0.,145.,RIN=100.,ROUT=5000.,ndots=15)
#profile_redMapper_indcat('gx_CS82_redMapper.fits','CS82_bin1',20.,23.42,RIN=100.,ROUT=5000.,ndots=15)
#profile_redMapper_indcat('gx_CS82_redMapper.fits','CS82_bin2',23.42,28.3,RIN=100.,ROUT=5000.,ndots=15)
#profile_redMapper_indcat('gx_CS82_redMapper.fits','CS82_bin3',28.3,39.7,RIN=100.,ROUT=5000.,ndots=15)


'''
print 'Now is trying to fit a NFW profile...'

try:
	nfw          = NFW_stack_fit(R[:BIN],shear[:BIN],err_et[:BIN],zmean,roc)
except:
	nfw          = [-999.,0.,-100.,[0.,0.],[0.,0.],-999.,0.]

c          = nfw[5]
CHI_nfw    = nfw[2]
R200       = nfw[0]
error_R200 = nfw[1]
x2         = nfw[3]
y2         = nfw[4]


fig = plt.figure(figsize=(8, 6))  #tamano del plot
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) #divide en 2 el eje x, en 1 el eje y y da la razon de alturas

#asigna los sublots

ax = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

#grafica


blancox=5000.
blancoy=5000.

name_label = 'redMapper'

x    = np.zeros(2,float)
y    = np.zeros(2,float)
x[1] = R.max()



ax.plot(R,shear,'ko')
ax.plot(blancox,blancoy,'w.',label=name_label)#

ax.legend(loc=1,frameon=False, scatterpoints = 1)
ax.plot(x2,y2,'C0-',label='NFW: c='+str('%.1f' % c)+', $R_{200}$ = '+str('%.2f' % R200)+' $\pm$ '+str('%.2f' % error_R200)+' Mpc$\,h^{-1}_{70}$, $\chi_{red}^{2} =$'+str('%.1f' % CHI_nfw)) #, c='+str(round(c)))
ax.errorbar(R, shear, yerr=err_et, fmt='None', ecolor='k')


#legend
matplotlib.rcParams['legend.fontsize'] = 11.
ax.legend(loc=1,frameon=False, scatterpoints = 1)


# axis detail
ax.axis([RIN/1000.,ROUT/1000.,1.,5000.])
ax.set_xscale('log')
ax.set_yscale('log')
ax.xaxis.set_ticks(np.arange(RIN/1000., ROUT/1000., 300.))
ax.set_xticklabels(np.arange(RIN/1000., ROUT/1000., 300.))
ax.yaxis.set_ticks(np.arange(10., 200., 100.))
ax.set_yticklabels(np.arange(10., 200., 100.))

#label					
ax.set_ylabel(u'$\Delta\Sigma_{\parallel} (M_{\odot}\,pc^{-2})$',fontsize=15)

#-----------------------------

ax2.plot(x,y,'k')
ax2.plot(R,cero,'kx')
ax2.errorbar(R,cero, yerr=err_ex, fmt='None', ecolor='k')

#axis details
ax2.axis([RIN/1000.,ROUT/1000.,-25.,30.])
ax2.yaxis.set_ticks(np.arange(-20., 20.,15.))
ax2.set_yticklabels(np.arange(-20., 20., 15.))
ax2.set_xscale('log')
ax2.xaxis.set_ticks(np.arange(RIN/1000.,ROUT/1000., 0.5))
ax2.set_xticklabels(np.arange(RIN/1000.,ROUT/1000., 0.5))



#labels
ax2.set_ylabel(r'$\Delta\Sigma_{\times} $',fontsize=15)
ax2.set_xlabel('r [$h^{-1}_{70}\,$Mpc]',fontsize=15)

#to join the plots
fig.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

plt.subplots_adjust(left=0.15, right=0.95, top=0.85, bottom=0.1)
'''
