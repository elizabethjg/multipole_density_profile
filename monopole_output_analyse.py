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
from profiles_fit import chi_red

# def plot_monopole_miss(folder,file_profile,out_file,m200_maria,pcc_maria):
def plot_monopole_miss(folder,file_profile,out_file):
	
	print file_profile
	
	file_mcmc = 'monopole_misscentred_'+file_profile
	
	os.system('mkdir '+folder+'plots_monopole_mcmc/')
	
	mcmc = (np.loadtxt(folder+file_mcmc)).T
	
	
	f, ax = plt.subplots(2, 1, figsize=(10,5))
	ax[0].plot(mcmc[0],'k.',alpha=0.3)
	ax[0].axvline(1500)
	# ax[0].axhline(np.log10(m200_maria*1.e14))
	ax[1].plot(mcmc[1],'k.',alpha=0.3)
	ax[1].axvline(1500)
	# ax[1].axhline(pcc_maria)
	plt.savefig(folder+'plots_monopole_mcmc/out_'+file_mcmc[:-3]+'png')
	
	mcmc = mcmc[:,1500:]
	
	labels = ['M200','pcc']
	
	
	
	m200 = np.percentile(mcmc[0], [16, 50, 84])
	pcc    = np.percentile(mcmc[1], [16, 50, 84])
	
	M200   = 10**(m200[1])/1.e14
	e_M200 = (10**(m200[1])*np.log(10.)*np.diff(m200))/1.e14
	
	print ' M200 ',M200,' +/- ',e_M200
	print ' pcc ',pcc[1],' +/- ',np.diff(pcc)
	
	
	fig = corner.corner(mcmc.T, labels=labels)
	plt.savefig(folder+'plots_monopole_mcmc/'+file_mcmc[:-3]+'png')
	
		
	profile = np.loadtxt(folder+file_profile).T
	
	f = open(folder+file_profile,'r')
	lines = f.readlines()
	j = lines[1].find('=')+1
	Mguess = (float(lines[1][j:-2])*1.e14)
	j = lines[2].find('=')+1
	zmean = float(lines[2][j:-2])
	
	r  = np.logspace(np.log10(min(profile[0])),
					np.log10(max(profile[0])),10)

	######################
	#COMPUTE CHI_SQUARE

	multipoles = multipole_shear_parallel(profile[0],M200=M200*1.e14,components = ['t'],
										misscentred = True,
										ellip=0.,z=zmean,
										verbose=True,ncores=4)
	
	model_t = model_Gamma(multipoles,'t', misscentred = True,pcc=pcc[1])

	chi_t = chi_red(model_t,profile[1],profile[2],1)

	#--------------------

	
	multipoles = multipole_shear_parallel(r,M200=M200*1.e14,components = ['t'],
										misscentred = True,
										ellip=0.,z=zmean,
										verbose=True,ncores=4)
	
	model_t = model_Gamma(multipoles,'t', misscentred = True,pcc=pcc[1])

	# multipoles = multipole_shear_parallel(r,M200=m200_maria*1.e14,components = ['t'],
										# misscentred = True,
										# ellip=0.,z=zmean,
										# verbose=True,ncores=4)
	
	# model_maria = model_Gamma(multipoles,'t', misscentred = True,pcc=pcc_maria)
	
	
	f, ax = plt.subplots(figsize=(6.5,5))
	ax.plot(profile[0],profile[1],'C0o')
	ax.plot(r,model_t,'C1',label = 'mcmc result')
	# ax.plot(r,model_maria,'C2',label = 'maria')
	ax.errorbar(profile[0],profile[1],yerr=profile[2],fmt = 'none',ecolor='C0')
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('R [mpc]')
	ax.set_ylim(1,300)
	ax.set_xlim(0.1,5.2)
	ax.xaxis.set_ticks([0.1,1,5])
	ax.set_xticklabels([0.1,1,5])
	ax.yaxis.set_ticks([1,10,100])
	ax.set_yticklabels([1,10,100])
	plt.legend()
	f.subplots_adjust(hspace=0,wspace=0)
	plt.savefig(folder+'plots_monopole_mcmc/'+'t'+file_profile[:-3]+'png')
	
	#/////// save file ///////
	
	f1=open(folder+out_file,'a')
	f1.write(str('%.2f' % (M200))+'   '+str('%.2f' % (e_M200[0]))+'   '+str('%.2f' % (e_M200[1]))+'   ')
	f1.write(str('%.2f' % (pcc[1]))+'   '+str('%.2f' % (np.diff(pcc)[0]))+'   '+str('%.2f' % (np.diff(pcc)[1]))+'   ')
	f1.write(str('%.2f' % chi_t)+'   \n')
	f1.close()
	
	f=open(folder+file_profile,'a')
	f.write('#'+str('%.2f' % (M200))+'   \n')
	f.write('#'+str('%.2f' % (pcc[1]))+'   \n')
	f.close()


folder        = u'/home/eli/Documentos/Astronomia/posdoc/halo-elongation/redMapper/profiles/'
# folder        = u'/home/eli/Documentos/PostDoc/halo-elongation/redMapper/profiles_original/'

os.system('rm '+folder+'out_table_monopole')
plot_monopole_miss(folder,'profile_total.cat','out_table_monopole')
plot_monopole_miss(folder,'profile_median_bin1.cat','out_table_monopole')
plot_monopole_miss(folder,'profile_median_bin2.cat','out_table_monopole')
# plot_monopole_miss(folder,'profile_terciles_bin1.cat','out_table_monopole')
# plot_monopole_miss(folder,'profile_terciles_bin2.cat','out_table_monopole')
# plot_monopole_miss(folder,'profile_terciles_bin3.cat','out_table_monopole')
plot_monopole_miss(folder,'profile_medianz1.cat','out_table_monopole')
plot_monopole_miss(folder,'profile_medianz2.cat','out_table_monopole')
# plot_monopole_miss(folder,'profile_terciles_z1.cat','out_table_monopole')
# plot_monopole_miss(folder,'profile_terciles_z2.cat','out_table_monopole')
# plot_monopole_miss(folder,'profile_terciles_z3.cat','out_table_monopole')

'''
maria_result = np.loadtxt(folder+'../../maria_result').T
m200_maria   = maria_result[0]/0.7
em200_maria  = maria_result[1]/0.7
pcc_maria    = maria_result[2]

plot_monopole_miss(folder,'profile_CS82_original_bin1.cat','out_CS82_table',m200_maria[0],pcc_maria[0])
plot_monopole_miss(folder,'profile_CS82_original_bin2.cat','out_CS82_table',m200_maria[1],pcc_maria[1])
plot_monopole_miss(folder,'profile_CS82_original_bin3.cat','out_CS82_table',m200_maria[2],pcc_maria[2])
plot_monopole_miss(folder,'profile_CS82_original_bin4.cat','out_CS82_table',m200_maria[3],pcc_maria[3])

lmean = [21.75,25.76,32.91,57.59]
'''
