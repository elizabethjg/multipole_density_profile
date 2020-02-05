import numpy as np
from matplotlib import *
from multipoles_shear import *

def make_plot(file_name,folder):
          
     f = open(folder+file_name,'r')
     lines = f.readlines()
     j = lines[1].find('=')+1
     M200 = float(lines[1][j:-2])*1.e14

     print '####################'
     print file_name
     print 'M200 ',M200/1.e14
     print '####################'

     
     profile = np.loadtxt(folder+file_name).T
     
     
     
     out = multipole_shear_parallel(profile[0],M200=M200,z=0.2,
                                   ellip=0.2,misscentred=False,
                                   ncores=2)
     
     Gt    = model_Gamma(out,'t',misscentred=False)     
     Gtcos1 = model_Gamma(out,'tcos',misscentred=False)     
     Gxsin1 = model_Gamma(out,'xsin',misscentred=False)     
     
     out = multipole_shear_parallel(profile[0],M200=M200,z=0.2,
                                   ellip=0.3,misscentred=False,
                                   ncores=2)
     
     
     Gtcos2 = model_Gamma(out,'tcos',misscentred=False)     
     Gxsin2 = model_Gamma(out,'xsin',misscentred=False)     
     
     
     f, ax = plt.subplots(1, 2, figsize=(10,5))
     ax[0].plot(profile[0],Gt,'C1--')
     ax[0].plot(profile[0],Gtcos1,'C1-')
     ax[0].plot(profile[0],Gtcos2,'C1-')
     ax[0].plot(profile[0],profile[1],'C0o')
     ax[0].errorbar(profile[0],profile[1],yerr=profile[2])
     # ax[0].axis([0,10,1,150.])
     ax[0].set_xscale('log')
     ax[0].set_yscale('log')
     ax[0].set_xlabel('R [mpc]')
     ax[0].set_ylim(1,200)
     ax[0].set_xlim(0.1,10)
     ax[0].xaxis.set_ticks([0.1,1,5])
     ax[0].set_xticklabels([0.1,1,5])
     ax[0].yaxis.set_ticks([1,10,100])
     ax[0].set_yticklabels([1,10,100])
     
     ax[1].plot(profile[0],Gxsin1,'C1-')
     ax[1].plot(profile[0],Gxsin2,'C1-')
     ax[1].plot(profile[0],profile[3],'C0o')
     ax[1].errorbar(profile[0],profile[3],yerr=profile[4])
     # ax[1].axis([0,10,-150,150.])
     ax[1].set_xscale('log')
     ax[1].set_xlabel('R [mpc]')
     ax[1].set_ylim(-200,200)
     ax[1].set_xlim(0.1,10)
     ax[1].xaxis.set_ticks([0.1,1,5])
     ax[1].set_xticklabels([0.1,1,5])
     ax[1].yaxis.set_ticks([-100,0,100])
     ax[1].set_yticklabels([-100,0,100])
     
     f.subplots_adjust(hspace=0,wspace=0)
     plt.savefig(folder+'plots/'+file_name[:-4]+'.png')

folder = '/home/eli/Documentos/PostDoc/halo-elongation/redMapper/profiles_terciles_z033/'
#folder = '/home/eli/Documentos/Astronomia/posdoc/halo-elongation/redMapper/profiles_terciles/'
     
f = open(folder+'list_names','r')
lines = f.readlines()

for line in lines:
     make_plot(line[:-1],folder)
