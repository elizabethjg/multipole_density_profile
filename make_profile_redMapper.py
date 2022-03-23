import sys
sys.path.append('/home/elizabeth/lens_codes_v3.7')
import numpy as np
# from pylab import *
from astropy.io import fits
#parameters
from astropy.cosmology import LambdaCDM
from maria_func import *
from make_profile import *
cosmo = LambdaCDM(H0=70., Om0=0.3, Ode0=0.7)
from multiprocessing import Pool
from multiprocessing import Process
#parameters
cvel = 299792458;   # Speed of light (m.s-1)
G    = 6.670e-11;   # Gravitational constant (m3.kg-1.s-2)
pc   = 3.085678e16; # 1 pc (m)
Msun = 1.989e30 # Solar mass (kg)

def profile_redMapper(sample,lmin,lmax,zmin = 0.1, zmax = 0.33,
                      z_back = 0.1, odds_min = 0.5,
                      RIN = 100., ROUT = 10000., ndots = 20.):
        
     try:
          lmin     = lmin.astype(float)
          lmax     = lmax.astype(float)
          zmin     = zmin.astype(float)
          zmax     = zmax.astype(float)
          z_back   = z_back.astype(float)
          odds_min = odds_min.astype(float)
          RIN      = RIN.astype(float)
          ROUT     = ROUT.astype(float)
          ndots    = int(ndots.astype(float))
     except:
          print 'not running in parallel'

     folder = '/mnt/clemente/lensing/redMaPPer/compressed/'
     
     cfht     = fits.open(folder+'gx_CFHT_redMapper.fits')[1].data
     kids     = fits.open(folder+'gx_KiDS_redMapper.fits')[1].data
     cs82     = fits.open(folder+'gx_CS82_redMapper.fits')[1].data
     rcsl     = fits.open(folder+'gx_RCSL_redMapper.fits')[1].data
     angles   = fits.open(folder+'angles_redMapper_forprofile.fits')[1].data
     clusters = fits.open(folder+'redmapper_dr8_public_v6.3_catalog.fits')[1].data
     borderid = np.loadtxt(folder+'redMapperID_border.list')
     
     mkids = (~np.in1d(kids.ID,cfht.ID))*(kids.Z_B < 0.9)
     mrcsl = (~np.in1d(rcsl.ID,cfht.ID))*(~np.in1d(rcsl.ID,cs82.ID))*(~np.in1d(rcsl.ID,kids.ID))*(rcsl.Z_B < 1.3)
     mcs82 = (~np.in1d(cs82.ID,cfht.ID))*(cs82.Z_B < 1.3)
     mcfht = (cfht.Z_B < 1.3)
     
     kids = kids[mkids]
     cs82 = cs82[mcs82]
     rcsl = rcsl[mrcsl]
     cfht = cfht[mcfht]
     
     # MATCH ANGLES  
     
     ID = np.concatenate((cfht.ID,kids.ID,cs82.ID,rcsl.ID))
     IDc = clusters.ID
     
     ides,index,c = np.unique(ID,return_index=True,return_counts=True)
     angles_in = angles[np.in1d(IDc,ID)]
     sindex = np.argsort(index)
     angles = np.repeat(angles_in[sindex],c[sindex])
     
     #-----------------------------------------------
     
     lamb = np.concatenate((cfht.LAMBDA,kids.LAMBDA,cs82.LAMBDA,rcsl.LAMBDA))
     
     
     Z_B  = np.concatenate((cfht.Z_B,kids.Z_B,cs82.Z_B,rcsl.Z_B))
     ODDS = np.concatenate((cfht.ODDS,kids.ODDS,cs82.ODDS,rcsl.ODDS))
     zlambda = np.concatenate((cfht.Z_LAMBDA,kids.Z_LAMBDA,cs82.Z_LAMBDA,rcsl.Z_LAMBDA))
     zspec   = np.concatenate((cfht.Z_SPEC,kids.Z_SPEC,cs82.Z_SPEC,rcsl.Z_SPEC))
     Z_c      = zspec
     Z_c[Z_c<0] = zlambda[Z_c<0]
     
     
     mask_back = (Z_B > (Z_c + z_back))*(ODDS >= odds_min)
     mask_lens = (Z_c >= zmin)*(Z_c < zmax)*(~np.in1d(ID,borderid))*(lamb >= lmin)*(lamb < lmax)
       
     mask = mask_back*mask_lens
     
     Nclusters = len(np.unique(ID[mask]))
     
     #del(mask_back)
     #del(mask_lens)
     #del(ODDS)
     #del(zlambda)
     #del(zspec)
     #del(Z_B)     
     #del(ID)
     

     # Cosmological distances

     Z_c    = Z_c[mask]
     lamb   = lamb[mask]              
     zmean    = (Z_c).mean()
     zdisp    = (Z_c).std()
     H        = cosmo.H(zmean).value/(1.0e3*pc) #H at z_pair s-1 
     roc      = (3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_pair (kg.m-3)
     roc_mpc  = roc*((pc*1.0e6)**3.0)
     D_ang    = cosmo.angular_diameter_distance(zmean)
     kpcscale = D_ang*(((1.0/3600.0)*np.pi)/180.0)*1000.0

     #del(Z_c)
     
     dls  = np.concatenate((cfht.DLS,kids.DLS,cs82.DLS,rcsl.DLS))[mask]
     ds   = np.concatenate((cfht.DS,kids.DS,cs82.DS,rcsl.DS))[mask]
     dl   = np.concatenate((cfht.DL,kids.DL,cs82.DL,rcsl.DL))[mask]
     
     
     KPCSCALE   = dl*(((1.0/3600.0)*np.pi)/180.0)*1000.0
     BETA_array = dls/ds
     beta       = BETA_array.mean()
     
     Dl = dl*1.e6*pc
     sigma_c = (((cvel**2.0)/(4.0*np.pi*G*Dl))*(1./BETA_array))*(pc**2/Msun)
     
     #del(dls)
     #del(dl)
     #del(ds)
     
     print 'BETA.mean',beta
     
     SIGMAC = (((cvel**2.0)/(4.0*np.pi*G*Dl.mean())))*(pc**2/Msun)
     
     print 'SIGMA_C', SIGMAC
     
     # Compute tangential and cross components
     
     ra     = np.concatenate((cfht.RAJ2000,kids.RAJ2000,cs82.RAJ2000,rcsl.RAJ2000))[mask]
     dec    = np.concatenate((cfht.DECJ2000,kids.DECJ2000,cs82.DECJ2000,rcsl.DECJ2000))[mask]
          
     ALFA0  = np.concatenate((cfht.RA,kids.RA,cs82.RA,rcsl.RA))[mask]
     DELTA0 = np.concatenate((cfht.DEC,kids.DEC,cs82.DEC,rcsl.DEC))[mask]
     
     rads, theta, test1,test2 = eq2p2(np.deg2rad(ra),
                              np.deg2rad(dec),
                              np.deg2rad(ALFA0),
                              np.deg2rad(DELTA0))
     
     #del(ra)
     #del(dec)
     #del(ALFA0)
     #del(DELTA0)
     
     
     theta2 = (2.*np.pi - theta) +np.pi/2.
     theta_ra = theta2
     theta_ra[theta2 > 2.*np.pi] = theta2[theta2 > 2.*np.pi] - 2.*np.pi

     #Correct polar angle for e1, e2
     theta = theta+np.pi/2.

     e1     = np.concatenate((cfht.e1,kids.e1,cs82.e1,rcsl.e1))[mask]
     e2     = np.concatenate((cfht.e2,kids.e2,cs82.e2,rcsl.e2))[mask]
     
     #get tangential ellipticities 
     et = (-e1*np.cos(2*theta)-e2*np.sin(2*theta))*sigma_c
     #get cross ellipticities
     ex = (-e1*np.sin(2*theta)+e2*np.cos(2*theta))*sigma_c
     
     #del(e1)
     #del(e2)
     
     r=np.rad2deg(rads)*3600*KPCSCALE
     #del(rads)
     
     peso = np.concatenate((cfht.weight,kids.weight,cs82.weight,rcsl.weight))[mask]
     peso = peso/(sigma_c**2) 
     m    = np.concatenate((cfht.m,kids.m,cs82.m,rcsl.m))[mask]
    
     
     
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
     
     f1=open(folder+'profile_'+sample+'.cat','w')
     f1.write('# Nclusters = '+str('%8i' % Nclusters)+' \n')
     f1.write('# M200 = '+str('%.2f' % (M200_NFW/1.e14))+' \n')
     f1.write('# z_mean = '+str('%.2f' % zmean)+' \n')
     f1.write('# z_back = '+str('%.2f' % z_back)+' \n')
     f1.write('# odds_min = '+str('%.1f' % odds_min)+' \n')
     f1.write('# l_min = '+str('%.1f' % lmin)+' \n')
     f1.write('# l_max = '+str('%.1f' % lmax)+' \n')
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
          
          f1=open(folder+'profile_'+sample+'.cat','w')
          f1.write('# Nclusters = '+str('%8i' % Nclusters)+' \n')
          f1.write('# M200 = '+str('%.2f' % (M200_NFW/1.e14))+' \n')
          f1.write('# z_mean = '+str('%.2f' % zmean)+' \n')
          f1.write('# z_back = '+str('%.2f' % z_back)+' \n')
          f1.write('# odds_min = '+str('%.1f' % odds_min)+' \n')
          f1.write('# l_min = '+str('%.1f' % lmin)+' \n')
          f1.write('# l_max = '+str('%.1f' % lmax)+' \n')
          f1.write('# l_mean = '+str('%.1f' % np.mean(lamb))+' \n')
          f1.write('# z_min = '+str('%.1f' % zmin)+' \n')
          f1.write('# z_max = '+str('%.1f' % zmax)+' \n')
          f1.write('# R,shear,err_et,cero,err_ex \n')
          profile = np.column_stack((profile[0]*1.0e-3,profile[1],profile[5],profile[2],profile[6]))
          np.savetxt(f1,profile,fmt = ['%12.6f']*5)
          f1.close()
     
     at     = theta_ra - angles.theta[mask]
     atwl   = theta_ra - angles.theta_wlum[mask]
     atwd   = theta_ra - angles.theta_wd[mask]
     atp    = theta_ra - angles.theta_pcut[mask]
     atpwl  = theta_ra - angles.theta_pcut_wlum[mask]
     atpwd  = theta_ra - angles.theta_pcut_wd[mask]
     atpwdl = theta_ra - angles.theta_pcut_wdl[mask]
     
     write_profile(at,sample+'_t')
     write_profile(atwl,sample+'_twl')
     write_profile(atwd,sample+'_twd')
     write_profile(atp,sample+'_tp')
     write_profile(atpwl,sample+'_tpwl')
     write_profile(atpwd,sample+'_tpwd')
     write_profile(atpwdl,sample+'_tpwdl')
     write_profile(theta_ra,sample+'_control')


def profile_redMapper_randomchoice(sample,lmin=20,lmax=150,zmin = 0.1, zmax = 0.4,
                      z_back = 0.1, odds_min = 0.5,
                      RIN = 100., ROUT = 5000., ndots = 10.):
        
     try:
          lmin     = lmin.astype(float)
          lmax     = lmax.astype(float)
          zmin     = zmin.astype(float)
          zmax     = zmax.astype(float)
          z_back   = z_back.astype(float)
          odds_min = odds_min.astype(float)
          RIN      = RIN.astype(float)
          ROUT     = ROUT.astype(float)
          ndots    = int(ndots.astype(float))
     except:
          print 'not running in parallel'

     folder = '/mnt/clemente/lensing/redMaPPer/'
     
     cfht     = fits.open(folder+'gx_CFHT_redMapper.fits')[1].data
     kids     = fits.open(folder+'gx_KiDS_redMapper.fits')[1].data
     cs82     = fits.open(folder+'gx_CS82_redMapper.fits')[1].data
     rcsl     = fits.open(folder+'gx_RCSL_redMapper.fits')[1].data
     angles   = fits.open(folder+'angles_redMapper_forprofile.fits')[1].data
     clusters = fits.open(folder+'redmapper_dr8_public_v6.3_catalog.fits')[1].data
     borderid = np.loadtxt(folder+'redMapperID_border.list')
     
     mkids = (~np.in1d(kids.ID,cfht.ID))*(kids.Z_B < 0.9)
     mrcsl = (~np.in1d(rcsl.ID,cfht.ID))*(~np.in1d(rcsl.ID,cs82.ID))*(~np.in1d(rcsl.ID,kids.ID))*(rcsl.Z_B < 1.3)
     mcs82 = (~np.in1d(cs82.ID,cfht.ID))*(cs82.Z_B < 1.3)
     mcfht = (cfht.Z_B < 1.3)
     
     kids = kids[mkids]
     cs82 = cs82[mcs82]
     rcsl = rcsl[mrcsl]
     cfht = cfht[mcfht]
     
     # MATCH ANGLES  
     
     ID = np.concatenate((cfht.ID,kids.ID,cs82.ID,rcsl.ID))
     IDc = clusters.ID
     
     ides,index,c = np.unique(ID,return_index=True,return_counts=True)
     angles_in = angles[np.in1d(IDc,ID)]
     sindex = np.argsort(index)
     angles = np.repeat(angles_in[sindex],c[sindex])
     
     #-----------------------------------------------
     
     lamb = np.concatenate((cfht.LAMBDA,kids.LAMBDA,cs82.LAMBDA,rcsl.LAMBDA))
     
     
     Z_B  = np.concatenate((cfht.Z_B,kids.Z_B,cs82.Z_B,rcsl.Z_B))
     ODDS = np.concatenate((cfht.ODDS,kids.ODDS,cs82.ODDS,rcsl.ODDS))
     zlambda = np.concatenate((cfht.Z_LAMBDA,kids.Z_LAMBDA,cs82.Z_LAMBDA,rcsl.Z_LAMBDA))
     zspec   = np.concatenate((cfht.Z_SPEC,kids.Z_SPEC,cs82.Z_SPEC,rcsl.Z_SPEC))
     Z_c      = zspec
     Z_c[Z_c<0] = zlambda[Z_c<0]
     
     
     mask_back = (Z_B > (Z_c + z_back))*(ODDS >= odds_min)
     mask_lens = (Z_c >= zmin)*(Z_c < zmax)*(~np.in1d(ID,borderid))*(lamb >= lmin)*(lamb < lmax)
     ides = np.unique(ID[mask_back*mask_lens])

     #del(mask_back)
     #del(mask_lens)
     del(ODDS)
     del(zlambda)
     del(zspec)
     del(Z_B)     
     Z_c0       = Z_c
     lamb0      = lamb

     for j in range(100):
     
          print 'muestra ',j
          t1 = time.time()
          
          id_choice = np.random.choice(ides,len(ides)/2,replace=False)
          mides     = np.in1d(ID,id_choice)
          
          mask = mask_back*mask_lens*mides
          
          Nclusters = len(np.unique(ID[mask]))
         
     
          # Cosmological distances
     
          Z_c    = Z_c0[mask]
          lamb   = lamb0[mask]              
          zmean    = (Z_c).mean()
          zdisp    = (Z_c).std()
          H        = cosmo.H(zmean).value/(1.0e3*pc) #H at z_pair s-1 
          roc      = (3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_pair (kg.m-3)
          roc_mpc  = roc*((pc*1.0e6)**3.0)
          D_ang    = cosmo.angular_diameter_distance(zmean)
          kpcscale = D_ang*(((1.0/3600.0)*np.pi)/180.0)*1000.0
     
          #del(Z_c)
          
          dls  = np.concatenate((cfht.DLS,kids.DLS,cs82.DLS,rcsl.DLS))[mask]
          ds   = np.concatenate((cfht.DS,kids.DS,cs82.DS,rcsl.DS))[mask]
          dl   = np.concatenate((cfht.DL,kids.DL,cs82.DL,rcsl.DL))[mask]
          
          
          KPCSCALE   = dl*(((1.0/3600.0)*np.pi)/180.0)*1000.0
          BETA_array = dls/ds
          beta       = BETA_array.mean()
          
          Dl = dl*1.e6*pc
          sigma_c = (((cvel**2.0)/(4.0*np.pi*G*Dl))*(1./BETA_array))*(pc**2/Msun)
          
          del(dls)
          del(dl)
          del(ds)
          
          print 'BETA.mean',beta
          
          SIGMAC = (((cvel**2.0)/(4.0*np.pi*G*Dl.mean())))*(pc**2/Msun)
          
          print 'SIGMA_C', SIGMAC
          
          # Compute tangential and cross components
          
          ra     = np.concatenate((cfht.RAJ2000,kids.RAJ2000,cs82.RAJ2000,rcsl.RAJ2000))[mask]
          dec    = np.concatenate((cfht.DECJ2000,kids.DECJ2000,cs82.DECJ2000,rcsl.DECJ2000))[mask]
               
          ALFA0  = np.concatenate((cfht.RA,kids.RA,cs82.RA,rcsl.RA))[mask]
          DELTA0 = np.concatenate((cfht.DEC,kids.DEC,cs82.DEC,rcsl.DEC))[mask]
          
          rads, theta, test1,test2 = eq2p2(np.deg2rad(ra),
                                   np.deg2rad(dec),
                                   np.deg2rad(ALFA0),
                                   np.deg2rad(DELTA0))
          
          del(ra)
          del(dec)
          del(ALFA0)
          del(DELTA0)
          
          
          theta2 = (2.*np.pi - theta) +np.pi/2.
          theta_ra = theta2
          theta_ra[theta2 > 2.*np.pi] = theta2[theta2 > 2.*np.pi] - 2.*np.pi
     
          #Correct polar angle for e1, e2
          theta = theta+np.pi/2.
     
          e1     = np.concatenate((cfht.e1,kids.e1,cs82.e1,rcsl.e1))[mask]
          e2     = np.concatenate((cfht.e2,kids.e2,cs82.e2,rcsl.e2))[mask]
          
          #get tangential ellipticities 
          et = (-e1*np.cos(2*theta)-e2*np.sin(2*theta))*sigma_c
          #get cross ellipticities
          ex = (-e1*np.sin(2*theta)+e2*np.cos(2*theta))*sigma_c
          
          del(e1)
          del(e2)
          
          r=np.rad2deg(rads)*3600*KPCSCALE
          del(rads)
          
          peso = np.concatenate((cfht.weight,kids.weight,cs82.weight,rcsl.weight))[mask]
          peso = peso/(sigma_c**2) 
          m    = np.concatenate((cfht.m,kids.m,cs82.m,rcsl.m))[mask]
     
          
          
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
          
          f1=open(folder+'test/profile_'+sample+'_'+str(int(j))+'.cat','w')
          f1.write('# Nclusters = '+str('%8i' % Nclusters)+' \n')
          f1.write('# M200 = '+str('%.2f' % (M200_NFW/1.e14))+' \n')
          f1.write('# z_mean = '+str('%.2f' % zmean)+' \n')
          f1.write('# z_back = '+str('%.2f' % z_back)+' \n')
          f1.write('# odds_min = '+str('%.1f' % odds_min)+' \n')
          f1.write('# l_min = '+str('%.1f' % lmin)+' \n')
          f1.write('# l_max = '+str('%.1f' % lmax)+' \n')
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
               
               f1=open(folder+'test/profile_'+sample+'_'+str(int(j))+'.cat','w')
               f1.write('# Nclusters = '+str('%8i' % Nclusters)+' \n')
               f1.write('# M200 = '+str('%.2f' % (M200_NFW/1.e14))+' \n')
               f1.write('# z_mean = '+str('%.2f' % zmean)+' \n')
               f1.write('# z_back = '+str('%.2f' % z_back)+' \n')
               f1.write('# odds_min = '+str('%.1f' % odds_min)+' \n')
               f1.write('# l_min = '+str('%.1f' % lmin)+' \n')
               f1.write('# l_max = '+str('%.1f' % lmax)+' \n')
               f1.write('# l_mean = '+str('%.1f' % np.mean(lamb))+' \n')
               f1.write('# z_min = '+str('%.1f' % zmin)+' \n')
               f1.write('# z_max = '+str('%.1f' % zmax)+' \n')
               f1.write('# R,shear,err_et,cero,err_ex \n')
               profile = np.column_stack((profile[0]*1.0e-3,profile[1],profile[5],profile[2],profile[6]))
               np.savetxt(f1,profile,fmt = ['%12.6f']*5)
               f1.close()
          
          atp    = theta_ra - angles.theta_pcut[mask]
          
          write_profile(atp,sample+'_tp')
          
          print (time.time()-t1)/60.


# profile_redMapper_randomchoice('medianas')

def profile_redMapper_indcat(survey,sample,lmin,lmax,zmin = 0.1, zmax = 0.33,
					z_back = 0.1, odds_min = 0.5,
					RIN = 100., ROUT = 10000., ndots = 20.,zlim = 1.3):
	
     print survey
     print sample
     
     sample = survey+'_'+sample
     
     try:
          lmin     = lmin.astype(float)
          lmax     = lmax.astype(float)
          zmin     = zmin.astype(float)
          zmax     = zmax.astype(float)
          z_back   = z_back.astype(float)
          odds_min = odds_min.astype(float)
          RIN      = RIN.astype(float)
          ROUT     = ROUT.astype(float)
          ndots    = int(ndots.astype(float))
     except:
          print 'not running in parallel'
     
     
     folder = '/mnt/clemente/lensing/redMaPPer/'
     
     backgx   = fits.open(folder+'gx_'+survey+'_redMapper.fits')[1].data
     angles   = fits.open(folder+'angles_redMapper_forprofile.fits')[1].data
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
     
     
     mask_back = (Z_B > (Z_c + z_back))*(ODDS >= odds_min)*(Z_B < zlim)
     mask_lens = (Z_c >= zmin)*(Z_c < zmax)*(~np.in1d(ID,borderid))*(lamb >= lmin)*(lamb < lmax)
     mask = mask_back*mask_lens
          
     mlambda = (lamb >= lmin)*(lamb < lmax) 
     
     mask = mask_back*mask_lens*mlambda
     
     Nclusters = len(np.unique(ID[mask]))
     
     del(mask_back)
     del(mask_lens)
     del(ODDS)
     del(zlambda)
     del(zspec)
     del(Z_B)     
     del(ID)
     
     
     # Cosmological distances
     
     Z_c    = Z_c[mask]
     lamb   = lamb[mask]              
     zmean    = (Z_c).mean()
     zdisp    = (Z_c).std()
     H        = cosmo.H(zmean).value/(1.0e3*pc) #H at z_pair s-1 
     roc      = (3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_pair (kg.m-3)
     roc_mpc  = roc*((pc*1.0e6)**3.0)
     D_ang    = cosmo.angular_diameter_distance(zmean)
     kpcscale = D_ang*(((1.0/3600.0)*np.pi)/180.0)*1000.0
     
     del(Z_c)
     
     dls  = backgx.DLS[mask]
     ds   = backgx.DS[mask]
     dl   = backgx.DL[mask]
     
     
     KPCSCALE   = dl*(((1.0/3600.0)*np.pi)/180.0)*1000.0
     BETA_array = dls/ds
     beta       = BETA_array.mean()
     
     Dl = dl*1.e6*pc
     sigma_c = (((cvel**2.0)/(4.0*np.pi*G*Dl))*(1./BETA_array))*(pc**2/Msun)
     
     del(dls)
     del(dl)
     del(ds)
     
     print 'BETA.mean',beta
     
     SIGMAC = (((cvel**2.0)/(4.0*np.pi*G*Dl.mean())))*(pc**2/Msun)
     
     print 'SIGMA_C', SIGMAC
     
     # Compute tangential and cross components
     
     ra     = backgx.RAJ2000[mask]
     dec    = backgx.DECJ2000[mask]
          
     ALFA0  = backgx.RA[mask]
     DELTA0 = backgx.DEC[mask]
     
     rads, theta, test1,test2 = eq2p2(np.deg2rad(ra),
                                   np.deg2rad(dec),
                                   np.deg2rad(ALFA0),
                                   np.deg2rad(DELTA0))
     
     del(ra)
     del(dec)
     del(ALFA0)
     del(DELTA0)
     
     
     theta2 = (2.*np.pi - theta) +np.pi/2.
     theta_ra = theta2
     theta_ra[theta2 > 2.*np.pi] = theta2[theta2 > 2.*np.pi] - 2.*np.pi
     
     #Correct polar angle for e1, e2
     theta = theta+np.pi/2.
     
     e1     = backgx.e1[mask]
     e2     = backgx.e2[mask]
     
     #get tangential ellipticities 
     et = (-e1*np.cos(2*theta)-e2*np.sin(2*theta))*sigma_c
     #get cross ellipticities
     ex = (-e1*np.sin(2*theta)+e2*np.cos(2*theta))*sigma_c
     
     del(e1)
     del(e2)
     
     r=np.rad2deg(rads)*3600*KPCSCALE
     del(rads)
     
     peso = backgx.weight[mask]
     peso = peso/(sigma_c**2) 
     m    = backgx.m[mask]
     
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
     
     f1=open(folder+'profile_'+sample+'.cat','w')
     f1.write('# Nclusters = '+str('%8i' % Nclusters)+' \n')
     f1.write('# M200 = '+str('%.2f' % (M200_NFW/1.e14))+' \n')
     f1.write('# z_mean = '+str('%.2f' % zmean)+' \n')
     f1.write('# z_back = '+str('%.2f' % z_back)+' \n')
     f1.write('# odds_min = '+str('%.1f' % odds_min)+' \n')
     f1.write('# l_min = '+str('%.1f' % lmin)+' \n')
     f1.write('# l_max = '+str('%.1f' % lmax)+' \n')
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
          
          f1=open(folder+'profile_'+sample+'.cat','w')
          f1.write('# Nclusters = '+str('%8i' % Nclusters)+' \n')
          f1.write('# M200 = '+str('%.2f' % (M200_NFW/1.e14))+' \n')
          f1.write('# z_mean = '+str('%.2f' % zmean)+' \n')
          f1.write('# z_back = '+str('%.2f' % z_back)+' \n')
          f1.write('# odds_min = '+str('%.1f' % odds_min)+' \n')
          f1.write('# l_min = '+str('%.1f' % lmin)+' \n')
          f1.write('# l_max = '+str('%.1f' % lmax)+' \n')
          f1.write('# l_mean = '+str('%.1f' % np.mean(lamb))+' \n')
          f1.write('# z_min = '+str('%.1f' % zmin)+' \n')
          f1.write('# z_max = '+str('%.1f' % zmax)+' \n')
          f1.write('# R,shear,err_et,cero,err_ex \n')
          profile = np.column_stack((profile[0]*1.0e-3,profile[1],profile[5],profile[2],profile[6]))
          np.savetxt(f1,profile,fmt = ['%12.6f']*5)
          f1.close()
     
     at     = theta_ra - angles.theta[mask]
     atwl   = theta_ra - angles.theta_wlum[mask]
     atwd   = theta_ra - angles.theta_wd[mask]
     atp    = theta_ra - angles.theta_pcut[mask]
     atpwl  = theta_ra - angles.theta_pcut_wlum[mask]
     atpwd  = theta_ra - angles.theta_pcut_wd[mask]
     atpwdl = theta_ra - angles.theta_pcut_wdl[mask]
     
     write_profile(at,sample+'_t')
     write_profile(atwl,sample+'_twl')
     write_profile(atwd,sample+'_twd')
     write_profile(atp,sample+'_tp')
     write_profile(atpwl,sample+'_tpwl')
     write_profile(atpwd,sample+'_tpwd')
     write_profile(atpwdl,sample+'_tpwdl')
     write_profile(theta_ra,sample+'_control')


def makeprofile_unpack(pinput):
	return profile_redMapper(*pinput)

def makeindprofile_unpack(pinput):
	return profile_redMapper_indcat(*pinput)

def makeprofile_parallel(samples,lmin,lmax,zmin,zmax,
                         z_back, odds_min,RIN, ROUT, ndots):
     	
     ncores = len(lmin)
          
     entrada = np.array([samples,lmin,lmax,zmin,zmax,
                         z_back,odds_min,RIN,ROUT,ndots]).T
                         
     pool = Pool(processes=(ncores))
     salida=np.array(pool.map(makeprofile_unpack, entrada))
     pool.terminate()

def makeindprofile_parallel(name_cat,samples,lmin,lmax,
                            zmin,zmax,z_back, odds_min,RIN, ROUT,
                            ndots,zlim):
     	
     ncores = len(lmin)
          
     entrada = np.array([name_cat,samples,lmin,lmax,zmin,zmax,
                         z_back,odds_min,RIN,ROUT,ndots,zlim]).T
                         
     pool = Pool(processes=(ncores))
     salida=np.array(pool.map(makeindprofile_unpack, entrada))
     pool.terminate()

#profile_redMapper('original_bin1',20.,23.42,RIN=100./0.7,ROUT=10000./0.7,ndots=20)
#profile_redMapper('original_bin2',23.42,28.3,RIN=100./0.7,ROUT=10000./0.7,ndots=20)
#profile_redMapper('original_bin3',28.3,39.7,RIN=100./0.7,ROUT=10000./0.7,ndots=20)
#profile_redMapper('original_bin4',39.7,145.,RIN=100./0.7,ROUT=10000./0.7,ndots=20)




