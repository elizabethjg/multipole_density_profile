import numpy as np
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 12})

def plot_profile(profile,label,color,ax1,ax2,ang=False):
     
     ax1.scatter(profile[0],profile[0]*profile[1],facecolor='none',edgecolors=color,label=label)
     ax1.errorbar(profile[0],profile[0]*profile[1],yerr=profile[0]*profile[2],color=color)
     ax1.set_xscale('log')
     ax1.axis([0,5,-30,30])
     if ang:
          ax1.text(0.2,-27.,ang)
     ax2.scatter(profile[0],profile[0]*profile[3],facecolor='none',edgecolors=color)
     ax2.errorbar(profile[0],profile[0]*profile[3],yerr=profile[0]*profile[4],color=color)
     ax2.set_xscale('log')
     ax2.axis([0,5,-30,30])
     
     f.subplots_adjust(hspace=0,wspace=0)


f, ax = plt.subplots(2,1, figsize=(4,5), sharex=True)
folder = u'/home/eli/Documentos/Astronomia/posdoc/halo-elongation/redMapper/profiles_ind/'

file_name    = 'profile_CFHT_total.cat'
profile_cfht = np.loadtxt(folder+file_name).T
file_name    = 'profile_CS82_total.cat'
profile_cs82 = np.loadtxt(folder+file_name).T
file_name    = 'profile_KiDS_total.cat'
profile_kids = np.loadtxt(folder+file_name).T
file_name    = 'profile_RCSL_total.cat'
profile_rcsl = np.loadtxt(folder+file_name).T

plot_profile(profile_cfht,'$CFHTLens$','C3',ax[0],ax[1])
plot_profile(profile_kids,'$KiDS-450$','C4',ax[0],ax[1])
plot_profile(profile_cs82,'$CS82$','C5',ax[0],ax[1])
plot_profile(profile_rcsl,'$RCSLens$','C6',ax[0],ax[1])
ax[0].axis([0,5,-10,60])
ax[0].legend()
ax[0].set_ylabel(r'$r \times \Delta\Sigma_{\parallel}$ ',fontsize = '14')
ax[1].set_ylabel(r'$r \times \Delta\Sigma_{\times}$    ',fontsize = '14')
ax[1].set_xlabel('$r [h^{-1}_{70}\,Mpc]$',fontsize = '14')
plt.rc('font', family='serif', size='12.0')
ax[1].xaxis.set_ticks([0.1,1.0,5.0])
ax[1].set_xticklabels([0.1,1.0,5.0])


plt.savefig(folder+'profile_compare_monopole.eps',format='eps',bbox_inches='tight')
###################################

f, ax = plt.subplots(4,4, figsize=(10,10), sharex=True, sharey=True)

file_name    = 'profile_CFHT_total_t.cat'
profile_cfht = np.loadtxt(folder+file_name).T
file_name    = 'profile_CS82_total_t.cat'
profile_cs82 = np.loadtxt(folder+file_name).T
file_name    = 'profile_KiDS_total_t.cat'
profile_kids = np.loadtxt(folder+file_name).T
file_name    = 'profile_RCSL_total_t.cat'
profile_rcsl = np.loadtxt(folder+file_name).T

plot_profile(profile_cfht,'$CFHTLens$','C3',ax[0,0],ax[1,0],ang = '$\phi_1$')
plot_profile(profile_kids,'$KiDS-450$','C4',ax[0,0],ax[1,0])
plot_profile(profile_cs82,'$CS82$','C5',ax[0,0],ax[1,0])
plot_profile(profile_rcsl,'$RCSLens$','C6',ax[0,0],ax[1,0])


file_name    = 'profile_CFHT_total_twl.cat'
profile_cfht = np.loadtxt(folder+file_name).T
file_name    = 'profile_CS82_total_twl.cat'
profile_cs82 = np.loadtxt(folder+file_name).T
file_name    = 'profile_KiDS_total_twl.cat'
profile_kids = np.loadtxt(folder+file_name).T
file_name    = 'profile_RCSL_total_twl.cat'
profile_rcsl = np.loadtxt(folder+file_name).T

plot_profile(profile_cfht,'CS82','C3',ax[0,1],ax[1,1],ang = '$\phi_L$')
plot_profile(profile_kids,'KiDS','C4',ax[0,1],ax[1,1])
plot_profile(profile_cs82,'KiDS','C5',ax[0,1],ax[1,1])
plot_profile(profile_rcsl,'KiDS','C6',ax[0,1],ax[1,1])

file_name    = 'profile_CFHT_total_twd.cat'
profile_cfht = np.loadtxt(folder+file_name).T
file_name    = 'profile_CS82_total_twd.cat'
profile_cs82 = np.loadtxt(folder+file_name).T
file_name    = 'profile_KiDS_total_twd.cat'
profile_kids = np.loadtxt(folder+file_name).T
file_name    = 'profile_RCSL_total_twd.cat'
profile_rcsl = np.loadtxt(folder+file_name).T

plot_profile(profile_cfht,'CS82','C3',ax[0,2],ax[1,2],ang = '$\phi_d$')
plot_profile(profile_kids,'KiDS','C4',ax[0,2],ax[1,2])
plot_profile(profile_cs82,'KiDS','C5',ax[0,2],ax[1,2])
plot_profile(profile_rcsl,'KiDS','C6',ax[0,2],ax[1,2])

file_name    = 'profile_CFHT_total_tp.cat'
profile_cfht = np.loadtxt(folder+file_name).T
file_name    = 'profile_CS82_total_tp.cat'
profile_cs82 = np.loadtxt(folder+file_name).T
file_name    = 'profile_KiDS_total_tp.cat'
profile_kids = np.loadtxt(folder+file_name).T
file_name    = 'profile_RCSL_total_tp.cat'
profile_rcsl = np.loadtxt(folder+file_name).T

empty_p = np.array([[],[],[],[],[]])
plot_profile(empty_p,'$CFHTLens$','C3',ax[0,3],ax[1,3])
plot_profile(empty_p,'$KiDS-450$','C4',ax[0,3],ax[1,3])
plot_profile(empty_p,'$CS82$','C5',ax[0,3],ax[1,3])
plot_profile(empty_p,'$RCSLens$','C6',ax[0,3],ax[1,3])
ax[0,3].legend()
ax[0,3].axis('off')
ax[1,3].axis('off')

plot_profile(profile_cfht,'CS82','C3',ax[2,0],ax[3,0],ang = '$\phi^*_1$')
plot_profile(profile_kids,'KiDS','C4',ax[2,0],ax[3,0])
plot_profile(profile_cs82,'KiDS','C5',ax[2,0],ax[3,0])
plot_profile(profile_rcsl,'KiDS','C6',ax[2,0],ax[3,0])

file_name    = 'profile_CFHT_total_tpwl.cat'
profile_cfht = np.loadtxt(folder+file_name).T
file_name    = 'profile_CS82_total_tpwl.cat'
profile_cs82 = np.loadtxt(folder+file_name).T
file_name    = 'profile_KiDS_total_tpwl.cat'
profile_kids = np.loadtxt(folder+file_name).T
file_name    = 'profile_RCSL_total_tpwl.cat'
profile_rcsl = np.loadtxt(folder+file_name).T

plot_profile(profile_cfht,'CS82','C3',ax[2,1],ax[3,1],ang = '$\phi^*_L$')
plot_profile(profile_kids,'KiDS','C4',ax[2,1],ax[3,1])
plot_profile(profile_cs82,'KiDS','C5',ax[2,1],ax[3,1])
plot_profile(profile_rcsl,'KiDS','C6',ax[2,1],ax[3,1])

file_name    = 'profile_CFHT_total_tpwd.cat'
profile_cfht = np.loadtxt(folder+file_name).T
file_name    = 'profile_CS82_total_tpwd.cat'
profile_cs82 = np.loadtxt(folder+file_name).T
file_name    = 'profile_KiDS_total_tpwd.cat'
profile_kids = np.loadtxt(folder+file_name).T
file_name    = 'profile_RCSL_total_tpwd.cat'
profile_rcsl = np.loadtxt(folder+file_name).T

plot_profile(profile_cfht,'CS82','C3',ax[2,2],ax[3,2],ang = '$\phi^*_d$')
plot_profile(profile_kids,'KiDS','C4',ax[2,2],ax[3,2])
plot_profile(profile_cs82,'KiDS','C5',ax[2,2],ax[3,2])
plot_profile(profile_rcsl,'KiDS','C6',ax[2,2],ax[3,2])

file_name    = 'profile_CFHT_total_control.cat'
profile_cfht = np.loadtxt(folder+file_name).T
file_name    = 'profile_CS82_total_control.cat'
profile_cs82 = np.loadtxt(folder+file_name).T
file_name    = 'profile_KiDS_total_control.cat'
profile_kids = np.loadtxt(folder+file_name).T
file_name    = 'profile_RCSL_total_control.cat'
profile_rcsl = np.loadtxt(folder+file_name).T

plot_profile(profile_cfht,'CS82','C3',ax[2,3],ax[3,3],ang = '$control$')
plot_profile(profile_kids,'KiDS','C4',ax[2,3],ax[3,3])
plot_profile(profile_cs82,'KiDS','C5',ax[2,3],ax[3,3])
plot_profile(profile_rcsl,'KiDS','C6',ax[2,3],ax[3,3])

ax[0,0].set_ylabel(r'$ r \times \Gamma_{t \cos{\theta}}$ ',fontsize = '14')
ax[1,0].set_ylabel(r'$ r \times \Gamma_{\times \sin{\theta}}$',fontsize = '14')
ax[2,0].set_ylabel(r'$ r \times \Gamma_{t \cos{\theta}}$ ',fontsize = '14')
ax[3,0].set_ylabel(r'$ r \times \Gamma_{\times \sin{\theta}}$',fontsize = '14')

ax[3,0].set_xlabel('$r [h^{-1}_{70}\,Mpc]$',fontsize='14')
ax[3,1].set_xlabel('$r [h^{-1}_{70}\,Mpc]$',fontsize='14')
ax[3,2].set_xlabel('$r [h^{-1}_{70}\,Mpc]$',fontsize='14')
ax[3,3].set_xlabel('$r [h^{-1}_{70}\,Mpc]$',fontsize='14')
plt.rc('font', family='serif', size='12.0')
ax[3,0].xaxis.set_ticks([0.2,1.0,4.0])
ax[3,1].xaxis.set_ticks([0.2,1.0,4.0])
ax[3,2].xaxis.set_ticks([0.2,1.0,4.0])
ax[3,3].xaxis.set_ticks([0.2,1.0,4.0])
ax[3,0].set_xticklabels([0.2,1.0,4.0])
ax[3,1].set_xticklabels([0.2,1.0,4.0])
ax[3,2].set_xticklabels([0.2,1.0,4.0])
ax[3,3].set_xticklabels([0.2,1.0,4.0])


plt.savefig(folder+'profile_compare_quadrupole.eps',format='eps',bbox_inches='tight')
###################################
