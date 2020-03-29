import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
sys.path.append('/home/elizabeth/multipole_density_profile')
sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
from quadrupole_output_analyse import *
import argparse
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

parser = argparse.ArgumentParser()
parser.add_argument('-folder', action='store', dest='folder',default='/mnt/clemente/lensing/redMaPPer/')
parser.add_argument('-file', action='store', dest='file_name', default='table_mcmc.out')
parser.add_argument('-profile', action='store', dest='profile', default='profile.cat')
parser.add_argument('-ncores', action='store', dest='ncores', default=20)
parser.add_argument('-misscentred', action='store', dest='miss', default='False')
args = parser.parse_args()

folder    = args.folder
out_file  = args.file_name
profile   = args.profile
ncores    = int(args.ncores)
if 'True' in args.miss:
	miss      = True
elif 'False' in args.miss:
	miss      = False

matplotlib.rcParams.update({'font.size': 14})

ft, axt = plt.subplots(2,6, figsize=(13,4.5), sharex=True)
# fx, axx = plt.subplots(2,3, figsize=(10,6), sharex=True, sharey=True)
ft.subplots_adjust(hspace=0,wspace=0)
# fx.subplots_adjust(hspace=0,wspace=0)

final_plots(folder,file_name,'t',False,2,axt[0,0],axt[1,0],r'$\phi_1$')
final_plots(folder,file_name,'twl',False,2,axt[0,1],axt[1,1],r'$\phi_L$')
final_plots(folder,file_name,'twd',False,2,axt[0,2],axt[1,2],r'$\phi_d$')
final_plots(folder,file_name,'tp',False,2,axt[0,3],axt[1,3],r'$\phi^*_1$')
final_plots(folder,file_name,'tpwl',False,2,axt[0,4],axt[1,4],r'$\phi^*_L$')
final_plots(folder,file_name,'tpwd',False,2,axt[0,5],axt[1,5],r'$\phi^*_d$')

plt.setp(axt[0,1].get_yticklabels(), visible=False)
plt.setp(axt[0,2].get_yticklabels(), visible=False)
plt.setp(axt[0,3].get_yticklabels(), visible=False)
plt.setp(axt[0,4].get_yticklabels(), visible=False)
plt.setp(axt[0,5].get_yticklabels(), visible=False)

plt.setp(axt[1,1].get_yticklabels(), visible=False)
plt.setp(axt[1,2].get_yticklabels(), visible=False)
plt.setp(axt[1,3].get_yticklabels(), visible=False)
plt.setp(axt[1,4].get_yticklabels(), visible=False)
plt.setp(axt[1,5].get_yticklabels(), visible=False)


plt.rc('font', family='serif', size='12.0')
axt[0,0].set_ylabel(r'$ \Gamma_{t \cos{\theta}} [h_{70}M_\odot/pc]$ ',labelpad=10)
axt[1,0].set_ylabel(r'$ \Gamma_{\times \sin{\theta}} [h_{70}M_\odot/pc]$',labelpad=10)
# axt[1,2].legend(loc=3,fontsize=11)
axt[1,5].legend(loc=4,fontsize=11)

# fx.savefig('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/redMapper/'+file_name[:-4]+'_xsin.pdf',bbox_inches='tight')
ft.savefig('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/redMapper/'+file_name[:-4]+'.pdf',bbox_inches='tight')

'''

plot_mcmc_quadrupole_out(folder,profile,'twl',miss,0,5000,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'twd',miss,0,5000,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'tp',miss,0,5000,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'tpwl',miss,0,5000,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'tpwd',miss,0,5000,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'control',miss,0,5000,out_file,ncores)

plot_mcmc_quadrupole_out(folder,profile,'t',miss,700,5000,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'twl',miss,700,5000,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'twd',miss,700,5000,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'tp',miss,700,5000,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'tpwl',miss,700,5000,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'tpwd',miss,700,5000,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'control',miss,700,5000,out_file,ncores)

plot_mcmc_quadrupole_out(folder,profile,'t',miss,0,700,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'twl',miss,0,700,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'twd',miss,0,700,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'tp',miss,0,700,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'tpwl',miss,0,700,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'tpwd',miss,0,700,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'control',miss,0,700,out_file,ncores)
'''
