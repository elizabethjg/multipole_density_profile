import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
sys.path.append('/home/elizabeth/multipole_density_profile')
sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
from quadrupole_output_analyse import *
import argparse

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

plot_mcmc_quadrupole_out_all(folder,file_name,out_file,ncores)

'''
plot_mcmc_quadrupole_out(folder,profile,'t',miss,0,5,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'twl',miss,0,5,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'twd',miss,0,5,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'tp',miss,0,5,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'tpwl',miss,0,5,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'tpwd',miss,0,5,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'control',miss,0,5,out_file,ncores)

plot_mcmc_quadrupole_out(folder,profile,'t',miss,1,5,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'twl',miss,1,5,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'twd',miss,1,5,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'tp',miss,1,5,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'tpwl',miss,1,5,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'tpwd',miss,1,5,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'control',miss,1,5,out_file,ncores)

plot_mcmc_quadrupole_out(folder,profile,'t',miss,0,1,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'twl',miss,0,1,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'twd',miss,0,1,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'tp',miss,0,1,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'tpwl',miss,0,1,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'tpwd',miss,0,1,out_file,ncores)
plot_mcmc_quadrupole_out(folder,profile,'control',miss,0,1,out_file,ncores)
'''
