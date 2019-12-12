import numpy as np
from matplotlib import *
from multipoles_shear import *

file_name = 'profile_bin4v2_tpwdl.cat'
folder = '/home/eli/Documentos/Astronomia/posdoc/halo-elongation/redMapper/'

profile = np.loadtxt(folder+file_name).T

f, ax = plt.subplots(1, 2, sharex=True, figsize=(10,5))

out = multipole_shear_parallel(profile[0],M200=3.e14,z=0.2,
                              ellip=0.2,misscentred=False,
                              ncores=2)

Gt1    = model_Gamma(out,'t',misscentred=False)     
Gtcos1 = model_Gamma(out,'tcos',misscentred=False)     
Gxsin1 = model_Gamma(out,'xsin',misscentred=False)     

out = multipole_shear_parallel(profile[0],M200=3.e14,z=0.2,
                              ellip=0.3,misscentred=False,
                              ncores=2)

Gt2    = model_Gamma(out,'t',misscentred=False)     
Gtcos2 = model_Gamma(out,'tcos',misscentred=False)     
Gxsin2 = model_Gamma(out,'xsin',misscentred=False)     



# ax[0].plot(profile[0],Gt,'C1-')
ax[0].plot(profile[0],Gtcos1,'C1-')
ax[0].plot(profile[0],Gtcos2,'C1-')
ax[0].plot(profile[0],profile[1],'C0o')
ax[0].errorbar(profile[0],profile[1],yerr=profile[2])
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_xlabel('R [mpc]')

# ax[1].plot(profile[0],Gxsin,'C1-')
ax[1].plot(profile[0],Gxsin1,'C1-')
ax[1].plot(profile[0],Gxsin2,'C1-')
ax[1].plot(profile[0],profile[3],'C0o')
ax[1].errorbar(profile[0],profile[3],yerr=profile[4])
ax[1].set_xscale('log')
ax[1].set_xlabel('R [mpc]')
f.subplots_adjust(hspace=0,wspace=0)
