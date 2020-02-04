import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
sys.path.append('/home/elizabeth/multipole_density_profile')
sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
import numpy as np
from pylab import *
from multipoles_shear import *
import time
import corner
import os

folder         = u'/home/eli/Documentos/PostDoc/halo-elongation/redMapper/profiles_terciles/'
file_profile   = 'profile_bin3_tpwl.cat'
out_file       = 'out_mcmc_table'

def plot_mcmc_quadrupole_out(folder,file_name,angle,miss,out_file,ncores):
	
	file_profile = file_name[:-4]+'_'+angle+'.cat'
	
	print file_profile

	if miss:
		file_mcmc_xsin = folder+'quadrupole_xsin_miss_'+file_name[:-4]+'_'+angle+'.out'
		file_mcmc_tcos = folder+'quadrupole_tcos_miss_'+file_name[:-4]+'_'+angle+'.out'
		plot_out       = file_name[:-4]+'_'+angle+'_miss.png'
	else:
		file_mcmc_xsin = folder+'quadrupole_xsin_'+file_name[:-4]+'_'+angle+'.out'
		file_mcmc_tcos = folder+'quadrupole_tcos_'+file_name[:-4]+'_'+angle+'.out'
		plot_out       = file_name[:-4]+'_'+angle+'.png'
	
	os.system('mkdir '+folder+'plots_mcmc/')
	
	f = open(folder+file_name,'r')
	lines = f.readlines()
	j = lines[2].find('=')+1
	zmean = float(lines[2][j:-2])
	pcc = float((lines[-1][1:-2]))
	M200 = float((lines[-2][1:-2]))*1.e14
	
	
	mcmc_xsin = (np.loadtxt(file_mcmc_xsin))
	mcmc_tcos = (np.loadtxt(file_mcmc_tcos))
	
	
	f, ax = plt.subplots(2, 1, figsize=(10,5))
	ax[0].plot(mcmc_xsin,'k.',alpha=0.3)
	ax[0].axvline(500)
	ax[1].plot(mcmc_tcos,'k.',alpha=0.3)
	ax[1].axvline(500)
	plt.savefig(folder+'plots_mcmc/mcmc_out_'+plot_out)
	
	mcmc_xsin = mcmc_xsin[:500] 
	mcmc_tcos = mcmc_tcos[:500] 
	
	e_x    = np.percentile(mcmc_xsin, [16, 50, 84])
	e_t    = np.percentile(mcmc_tcos, [16, 50, 84])
	
	#/////// xsin ///////
	f, ax = plt.subplots(1, 2, figsize=(10,5))
	
	
	ax[0].hist(mcmc_tcos,histtype='step')
	ax[0].set_xlabel('e_t')
	ax[1].hist(mcmc_xsin,histtype='step')
	ax[1].set_xlabel('e_x')
	plt.savefig(folder+'plots_mcmc/hist_'+plot_out)
	
	#/////// save file ///////
	
	f1=open(folder+out_file,'a')
	f1.write(file_profile+'  '+'  tcos  ') 
	f1.write(str('%.2f' % (e_t[1]))+'   '+str('%.2f' % (np.diff(e_t)[0]))+'   '+str('%.2f' % (np.diff(e_t)[1]))+'   ')
	f1.write('  xsin  ') 
	f1.write(str('%.2f' % (e_x[1]))+'   '+str('%.2f' % (np.diff(e_x)[0]))+'   '+str('%.2f' % (np.diff(e_x)[1]))+'   ')
	f1.close()
	
	profile = np.loadtxt(folder+file_profile).T
	
	f = open(folder+file_profile,'r')
	lines = f.readlines()
	j = lines[1].find('=')+1
	Mguess = (float(lines[1][j:-2])*1.e14)
	j = lines[2].find('=')+1
	zmean = float(lines[2][j:-2])
	
	r  = np.logspace(np.log10(min(profile[0])),
					np.log10(max(profile[0])),10)
	
	multipoles = multipole_shear_parallel(r,M200=M200,misscentred = miss,
								ellip=e_t[1],z=zmean,components = ['tcos','xsin'],
								verbose=False,ncores=ncores)
	
	modelt_t = model_Gamma(multipoles,'tcos', misscentred = miss, pcc = pcc)
	modelt_x = model_Gamma(multipoles,'xsin', misscentred = miss, pcc = pcc)
	
	multipoles = multipole_shear_parallel(r,M200=M200,misscentred = miss,
								ellip=e_x[1],z=zmean,components = ['tcos','xsin'],
								verbose=False,ncores=ncores)
	
	modelx_t = model_Gamma(multipoles,'tcos', misscentred = False)
	modelx_x = model_Gamma(multipoles,'xsin', misscentred = False)
	
	
	f, ax = plt.subplots(1, 2, figsize=(6.5,5))
	ax[0].plot(profile[0],profile[1],'k')
	ax[0].plot(r,modelt_t,'C2',label = 'tcos')
	ax[0].plot(r,modelx_t,'C3',label = 'xsin')
	ax[0].errorbar(profile[0],profile[1],yerr=profile[2],fmt = 'none',ecolor='k')
	ax[0].set_xscale('log')
	ax[0].set_yscale('log')
	ax[0].set_xlabel('R [mpc]')
	ax[0].set_ylim(0.1,200)
	ax[0].set_xlim(0.1,5.2)
	ax[0].xaxis.set_ticks([0.1,1,5])
	ax[0].set_xticklabels([0.1,1,5])
	ax[0].yaxis.set_ticks([1,10,100])
	ax[0].set_yticklabels([1,10,100])
	plt.legend()
	
	ax[1].plot(profile[0],profile[3],'k')
	ax[1].plot(r,modelt_x,'C2',label = 'tcos')
	ax[1].plot(r,modelx_x,'C3',label = 'xsin')
	ax[1].errorbar(profile[0],profile[3],yerr=profile[4],fmt = 'none',ecolor='k')
	ax[1].set_xlabel('R [mpc]')
	ax[1].set_xscale('log')
	ax[1].set_ylim(-150,150)
	ax[1].set_xlim(0.1,5.2)
	ax[1].xaxis.set_ticks([0.1,1,5])
	ax[1].set_xticklabels([0.1,1,5])
	ax[1].yaxis.set_ticks([-100,0,100])
	ax[1].set_yticklabels([-100,0,100])
	plt.legend()
	plt.savefig(folder+'plots_mcmc/'+plot_out)


folder = u'/home/eli/Documentos/PostDoc/halo-elongation/redMapper/profiles_terciles/'

os.system('rm '+folder+'out_mcmc_centred')
plot_mcmc_quadrupole_out(folder,'profile_terciles_total.cat','t',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_total.cat','twl',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_total.cat','twd',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_total.cat','tp',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_total.cat','tpwl',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_total.cat','tpwd',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_total.cat','tpwdl',False,'out_mcmc_centred',4)

plot_mcmc_quadrupole_out(folder,'profile_terciles_bin1.cat','t',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_bin1.cat','twl',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_bin1.cat','twd',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_bin1.cat','tp',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_bin1.cat','tpwl',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_bin1.cat','tpwd',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_bin1.cat','tpwdl',False,'out_mcmc_centred',4)

plot_mcmc_quadrupole_out(folder,'profile_terciles_bin2.cat','t',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_bin2.cat','twl',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_bin2.cat','twd',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_bin2.cat','tp',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_bin2.cat','tpwl',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_bin2.cat','tpwd',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_bin2.cat','tpwdl',False,'out_mcmc_centred',4)

plot_mcmc_quadrupole_out(folder,'profile_terciles_bin3.cat','t',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_bin3.cat','twl',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_bin3.cat','twd',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_bin3.cat','tp',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_bin3.cat','tpwl',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_bin3.cat','tpwd',False,'out_mcmc_centred',4)
plot_mcmc_quadrupole_out(folder,'profile_terciles_bin3.cat','tpwdl',False,'out_mcmc_centred',4)
