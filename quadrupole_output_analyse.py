import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
sys.path.append('/home/elizabeth/multipole_density_profile')
sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
import numpy as np
import matplotlib
# matplotlib.use('Agg')
from pylab import *
from multipoles_shear import *
import time
import os
from profiles_fit import chi_red

def plot_mcmc_quadrupole_out(folder,file_name,ang,miss,
							 rin,rout,out_file,ncores):
	
	file_profile = file_name[:-4]+'_'+ang+'.cat'
	print rin, rout
	print file_profile

	if miss:
		file_mcmc_xsin = folder+'quadrupole_xsin_miss_'+file_name[:-4]+'_'+ang+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		file_mcmc_tcos = folder+'quadrupole_tcos_miss_'+file_name[:-4]+'_'+ang+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		file_mcmc_both = folder+'quadrupole_both_miss_'+file_name[:-4]+'_'+ang+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		plot_out       = file_name[:-4]+'_'+ang+'_'+str(int(rin))+'_'+str(int(rout))+'_miss.png'
	else:
		file_mcmc_xsin = folder+'quadrupole_xsin_'+file_name[:-4]+'_'+ang+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		file_mcmc_tcos = folder+'quadrupole_tcos_'+file_name[:-4]+'_'+ang+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		file_mcmc_both = folder+'quadrupole_both_'+file_name[:-4]+'_'+ang+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		plot_out       = file_name[:-4]+'_'+ang+'_'+str(int(rin))+'_'+str(int(rout))+'.png'
	
	os.system('mkdir '+folder+'plots_mcmc_centred/')
	
	f = open(folder+file_name,'r')
	lines = f.readlines()
	j = lines[2].find('=')+1
	zmean = float(lines[2][j:-2])
	pcc = float((lines[-1][1:-2]))
	M200 = float((lines[-2][1:-2]))*1.e14
	
	try:
		mcmc_xsin = (np.loadtxt(file_mcmc_xsin))
		mcmc_tcos = (np.loadtxt(file_mcmc_tcos))
		mcmc_both = (np.loadtxt(file_mcmc_both))
	except:
		print file_mcmc_xsin
		print file_mcmc_tcos
		print file_mcmc_both
		f1=open(folder+out_file,'a')
		f1.write(file_profile+' ') 
		f1.write(str('%.2f' % (M200/1.e14))+'   '+str('%.2f' % (zmean))+'   '+str(int(rin))+'   '+str(int(rout))+'    ')
		f1.write('0.     0.01     0.01    0.     0.01     0.01    0.     0.01     0.01 \n')
		f1.close()
		print 'FILE NOT FOUND'
		# return None
	
	
	f, ax = plt.subplots(3, 1, figsize=(6,5))
	ax[0].plot(mcmc_tcos,'k.',alpha=0.3)
	ax[0].axvline(500)
	ax[1].plot(mcmc_xsin,'k.',alpha=0.3)
	ax[1].axvline(500)
	ax[2].plot(mcmc_both,'k.',alpha=0.3)
	ax[2].axvline(500)
	f.subplots_adjust(hspace=0,wspace=0)
	plt.savefig(folder+'plots_mcmc_centred/mcmc_out_'+plot_out)
	
	mcmc_xsin = mcmc_xsin[500:] 
	mcmc_tcos = mcmc_tcos[500:] 
	mcmc_both = mcmc_both[500:] 
	
	e_x    = np.percentile(mcmc_xsin, [16, 50, 84])
	e_t    = np.percentile(mcmc_tcos, [16, 50, 84])
	e_b    = np.percentile(mcmc_both, [16, 50, 84])
	
	f, ax = plt.subplots(figsize=(3.5,5))
	
	ax.hist(mcmc_tcos,histtype='step',label='tcos')
	ax.hist(mcmc_xsin,histtype='step',label='xsin')
	ax.hist(mcmc_both,histtype='step',label='both')
	ax.set_xlabel(r'$\epsilon$')
	f.subplots_adjust(hspace=0,wspace=0)
	plt.legend()
	plt.savefig(folder+'plots_mcmc_centred/hist_'+plot_out)
	
	
	profile = np.loadtxt(folder+file_profile).T
	
	f = open(folder+file_profile,'r')
	lines = f.readlines()
	j = lines[1].find('=')+1
	Mguess = (float(lines[1][j:-2])*1.e14)
	j = lines[2].find('=')+1
	zmean = float(lines[2][j:-2])

	######################
	#COMPUTE CHI_SQUARE

	multipoles = multipole_shear_parallel(profile[0],M200=M200,misscentred = miss,
								ellip=e_t[1],z=zmean,components = ['tcos'],
								verbose=False,ncores=ncores)
	
	modelt_t = model_Gamma(multipoles,'tcos', misscentred = miss, pcc = pcc)

	chi_t = chi_red(modelt_t,profile[1],profile[2],1)

	#--------------------

	multipoles = multipole_shear_parallel(profile[0],M200=M200,misscentred = miss,
								ellip=e_t[1],z=zmean,components = ['xsin'],
								verbose=False,ncores=ncores)
	
	modelx_x = model_Gamma(multipoles,'xsin', misscentred = miss, pcc = pcc)

	chi_x = chi_red(modelx_x,profile[3],profile[4],1)

	#--------------------

	multipoles = multipole_shear_parallel(profile[0],M200=M200,misscentred = miss,
								ellip=e_b[1],z=zmean,components = ['tcos','xsin'],
								verbose=False,ncores=ncores)
	
	modelb_t = model_Gamma(multipoles,'tcos', misscentred = miss, pcc = pcc)
	modelb_x = model_Gamma(multipoles,'xsin', misscentred = miss, pcc = pcc)

	chi_b = chi_red(np.append(modelb_t,modelb_x),
			np.append(profile[1],profile[3]),
			np.append(profile[2],profile[4]),1)

	######################


	
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
	
	modelx_t = model_Gamma(multipoles,'tcos', misscentred = miss, pcc = pcc)
	modelx_x = model_Gamma(multipoles,'xsin', misscentred = miss, pcc = pcc)

	multipoles = multipole_shear_parallel(r,M200=M200,misscentred = miss,
								ellip=e_b[1],z=zmean,components = ['tcos','xsin'],
								verbose=False,ncores=ncores)
	
	modelb_t = model_Gamma(multipoles,'tcos', misscentred = miss, pcc = pcc)
	modelb_x = model_Gamma(multipoles,'xsin', misscentred = miss, pcc = pcc)
	
	
	f, ax = plt.subplots(1, 2, figsize=(8,5))
	ax[0].scatter(profile[0],profile[1],facecolor='none',edgecolors='0.4')
	ax[0].plot(r,modelt_t,'C2',label = 'tcos')
	ax[0].fill_between(r, modelt_t - modelt_t*np.diff(e_t)[0], modelt_t + modelt_t*np.diff(e_t)[1],facecolor = 'C2', alpha = 0.3)
	ax[0].plot(r,modelx_t,'C3',label = 'xsin')
	ax[0].fill_between(r, modelx_t - modelx_t*np.diff(e_x)[0], modelx_t + modelx_t*np.diff(e_x)[1],facecolor = 'C3', alpha = 0.3)
	ax[0].plot(r,modelb_t,'C4',label = 'both')
	ax[0].fill_between(r, modelb_t - modelb_t*np.diff(e_b)[0], modelb_t + modelb_t*np.diff(e_b)[1],facecolor = 'C4', alpha = 0.3)
	ax[0].errorbar(profile[0],profile[1],yerr=profile[2],fmt = 'none',ecolor='0.4')
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
	
	ax[1].scatter(profile[0],profile[3],facecolor='none',edgecolors='0.4')
	ax[1].plot(r,modelt_x,'C2',label = 'tcos')
	ax[1].fill_between(r, modelt_x - modelt_x*np.diff(e_t)[0], modelt_x + modelt_x*np.diff(e_t)[1],facecolor = 'C2', alpha = 0.3)
	ax[1].plot(r,modelx_x,'C3',label = 'xsin')
	ax[1].fill_between(r, modelx_x - modelx_x*np.diff(e_x)[0], modelx_x + modelx_x*np.diff(e_x)[1],facecolor = 'C3', alpha = 0.3)
	ax[1].plot(r,modelb_x,'C4',label = 'both')
	ax[1].fill_between(r, modelb_x - modelb_x*np.diff(e_b)[0], modelb_x + modelb_x*np.diff(e_b)[1],facecolor = 'C4', alpha = 0.3)
	ax[1].errorbar(profile[0],profile[3],yerr=profile[4],fmt = 'none',ecolor='0.4')
	ax[1].set_xlabel('R [mpc]')
	ax[1].set_xscale('log')
	ax[1].set_ylim(-50,50)
	ax[1].set_xlim(0.1,5.2)
	ax[1].xaxis.set_ticks([0.1,1,5])
	ax[1].set_xticklabels([0.1,1,5])
	ax[1].yaxis.set_ticks([-30,-15,0,15,30])
	ax[1].set_yticklabels([-30,-15,0,15,30])
	plt.legend()
	f.subplots_adjust(hspace=0,wspace=0)
	plt.savefig(folder+'plots_mcmc_centred/'+plot_out)
	
	# #/////// save file ///////
	
	f1=open(folder+out_file,'a')
	f1.write(file_profile+' ') 
	f1.write(str('%.2f' % (M200/1.e14))+'   '+str('%.2f' % (zmean))+'   '+str(int(rin))+'   '+str(int(rout))+'    ')
	f1.write(str('%.2f' % (e_t[1]))+'   '+str('%.2f' % (np.diff(e_t)[0]))+'   '+str('%.2f' % (np.diff(e_t)[1]))+'   '+str('%.2f' % (chi_t))+'   ')
	f1.write(str('%.2f' % (e_x[1]))+'   '+str('%.2f' % (np.diff(e_x)[0]))+'   '+str('%.2f' % (np.diff(e_x)[1]))+'   '+str('%.2f' % (chi_x))+'   ')
	f1.write(str('%.2f' % (e_b[1]))+'   '+str('%.2f' % (np.diff(e_b)[0]))+'   '+str('%.2f' % (np.diff(e_b)[1]))+'   '+str('%.2f' % (chi_b))+' \n')
	f1.close()

def plot_mcmc_quadrupole_out_onlyboth(folder,file_name,ang,miss,
							 rin,rout,out_file,ncores):
	
	file_profile = file_name[:-4]+'_'+ang+'.cat'
	
	print file_profile

	if miss:
		file_mcmc_xsin = folder+'quadrupole_xsin_miss_'+file_name[:-4]+'_'+ang+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		file_mcmc_tcos = folder+'quadrupole_tcos_miss_'+file_name[:-4]+'_'+ang+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		file_mcmc_both = folder+'quadrupole_both_miss_'+file_name[:-4]+'_'+ang+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		plot_out       = file_name[:-4]+'_'+ang+'_'+str(int(rin))+'_'+str(int(rout))+'_miss.png'
	else:
		file_mcmc_xsin = folder+'quadrupole_xsin_'+file_name[:-4]+'_'+ang+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		file_mcmc_tcos = folder+'quadrupole_tcos_'+file_name[:-4]+'_'+ang+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		file_mcmc_both = folder+'quadrupole_both_'+file_name[:-4]+'_'+ang+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		plot_out       = file_name[:-4]+'_'+ang+'_'+str(int(rin))+'_'+str(int(rout))+'.png'
	
	os.system('mkdir '+folder+'plots_mcmc/')
	
	f = open(folder+file_name,'r')
	lines = f.readlines()
	j = lines[2].find('=')+1
	zmean = float(lines[2][j:-2])
	pcc = float((lines[-1][1:-2]))
	M200 = float((lines[-2][1:-2]))*1.e14
	
	try:
		mcmc_both = (np.loadtxt(file_mcmc_both))
	except:
		f1=open(folder+out_file,'a')
		f1.write(file_profile+' ') 
		f1.write(str('%.2f' % (M200/1.e14))+'   '+str('%.2f' % (zmean))+'   '+str(int(rin))+'   '+str(int(rout))+'    ')
		f1.write('0.     0.01     0.01 \n')
		f1.close()
		print 'FILE NOT FOUND'
		return None
	
	
	f, ax = plt.subplots(figsize=(5,2))
	ax.plot(mcmc_both,'k.',alpha=0.3)
	ax.axvline(500)
	f.subplots_adjust(hspace=0,wspace=0)
	plt.savefig(folder+'plots_mcmc/mcmc_out_'+plot_out)
	
	mcmc_both = mcmc_both[500:1000] 
	
	e_b    = np.percentile(mcmc_both, [16, 50, 84])
	
	f, ax = plt.subplots(figsize=(3,5))
	ax.hist(mcmc_both,histtype='step')
	ax.set_xlabel('e_b')
	f.subplots_adjust(hspace=0,wspace=0)
	plt.savefig(folder+'plots_mcmc/hist_'+plot_out)
	
	
	profile = np.loadtxt(folder+file_profile).T
	
	f = open(folder+file_profile,'r')
	lines = f.readlines()
	j = lines[1].find('=')+1
	Mguess = (float(lines[1][j:-2])*1.e14)
	j = lines[2].find('=')+1
	zmean = float(lines[2][j:-2])

	######################
	#COMPUTE CHI_SQUARE

	multipoles = multipole_shear_parallel(profile[0],M200=M200,misscentred = miss,
								ellip=e_b[1],z=zmean,components = ['tcos','xsin'],
								verbose=False,ncores=ncores)
	
	modelb_t = model_Gamma(multipoles,'tcos', misscentred = miss, pcc = pcc)
	modelb_x = model_Gamma(multipoles,'xsin', misscentred = miss, pcc = pcc)

	chi = chi_red(np.append(modelb_t,modelb_x),
			np.append(profile[1],profile[3]),
			np.append(profile[2],profile[4]),1)

	######################
	
	r  = np.logspace(np.log10(min(profile[0])),
					np.log10(max(profile[0])),15)
	

	multipoles = multipole_shear_parallel(r,M200=M200,misscentred = miss,
								ellip=e_b[1],z=zmean,components = ['tcos','xsin'],
								verbose=False,ncores=ncores)
	
	modelb_t = model_Gamma(multipoles,'tcos', misscentred = miss, pcc = pcc)
	modelb_x = model_Gamma(multipoles,'xsin', misscentred = miss, pcc = pcc)
	
	
	f, ax = plt.subplots(1, 2, figsize=(8,5))
	ax[0].scatter(profile[0],profile[1],facecolor='none',edgecolors='0.4')
	ax[0].plot(r,modelb_t,'C4')
	ax[0].errorbar(profile[0],profile[1],yerr=profile[2],fmt = 'none',ecolor='0.4')
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
	
	ax[1].scatter(profile[0],profile[3],facecolor='none',edgecolors='0.4')
	ax[1].plot(r,modelb_x,'C4')
	ax[1].errorbar(profile[0],profile[3],yerr=profile[4],fmt = 'none',ecolor='0.4')
	ax[1].set_xlabel('R [mpc]')
	ax[1].set_xscale('log')
	ax[1].set_ylim(-50,50)
	ax[1].set_xlim(0.1,5.2)
	ax[1].xaxis.set_ticks([0.1,1,5])
	ax[1].set_xticklabels([0.1,1,5])
	ax[1].yaxis.set_ticks([-30,-15,0,15,30])
	ax[1].set_yticklabels([-30,-15,0,15,30])
	plt.legend()
	f.subplots_adjust(hspace=0,wspace=0)
	plt.savefig(folder+'plots_mcmc/'+plot_out)
	
	#/////// save file ///////
	
	f1=open(folder+out_file,'a')
	f1.write(file_profile+' ') 
	f1.write(str('%.2f' % (M200/1.e14))+'   '+str('%.2f' % (zmean))+'   '+str(int(rin))+'   '+str(int(rout))+'    ')
	f1.write(str('%.2f' % (e_b[1]))+'   '+str('%.2f' % (np.diff(e_b)[0]))+'   '+str('%.2f' % (np.diff(e_b)[1]))+'   ')
	f1.write(str('%.2f' % chi)+'\n')
	f1.close()

def plot_mcmc_quadrupole_out_all(folder,file_name,out_file,ncores=20):
	
		
	
	angs = ['t','twl','twd','tp','tpwl',
              'tpwd','control']
	
	for ang in angs:
		
		t1 = time.time()
		
		file_profile = file_name[:-4]+'_'+ang+'.cat'
		
		print file_profile

		file_mcmc_both_miss = folder+'quadrupole_both_miss_'+file_name[:-4]+'_'+ang+'_0_5.out'
		file_mcmc_both = folder+'quadrupole_both_'+file_name[:-4]+'_'+ang+'_0_5.out'
		file_mcmc_both_miss_in = folder+'quadrupole_both_miss_'+file_name[:-4]+'_'+ang+'_0_1.out'
		file_mcmc_both_in = folder+'quadrupole_both_'+file_name[:-4]+'_'+ang+'_0_1.out'
		file_mcmc_both_miss_out = folder+'quadrupole_both_miss_'+file_name[:-4]+'_'+ang+'_1_5.out'
		file_mcmc_both_out = folder+'quadrupole_both_'+file_name[:-4]+'_'+ang+'_1_5.out'
		plot_out       = file_name[:-4]+'_'+ang+'.png'
		
		os.system('mkdir '+folder+'plots_all_mcmc/')
		
		f = open(folder+file_name,'r')
		lines = f.readlines()
		j = lines[2].find('=')+1
		zmean = float(lines[2][j:-2])
		pcc = float((lines[-1][1:-2]))
		M200 = float((lines[-2][1:-2]))*1.e14
		
		try:
			mcmc_both          = (np.loadtxt(file_mcmc_both))
			mcmc_both_miss     = (np.loadtxt(file_mcmc_both_miss))
			mcmc_both_in       = (np.loadtxt(file_mcmc_both_in))
			mcmc_both_miss_in  = (np.loadtxt(file_mcmc_both_miss_in))
			mcmc_both_out      = (np.loadtxt(file_mcmc_both_out))
			mcmc_both_miss_out = (np.loadtxt(file_mcmc_both_miss_out))
		except:
			
			f1=open(folder+out_file,'a')
			f1.write(file_profile+' ') 
			f1.write(str('%.2f' % (M200/1.e14))+'   '+str('%.2f' % (zmean))+'   0   5    ')
			f1.write('0.     0.01     0.01  0.0  0.     0.01     0.01  0.0 \n')
	
			f1.write(file_profile+' ') 
			f1.write(str('%.2f' % (M200/1.e14))+'   '+str('%.2f' % (zmean))+'   0   1    ')
			f1.write('0.     0.01     0.01  0.0  0.     0.01     0.01  0.0 \n')
	
			f1.write(file_profile+' ') 
			f1.write('0.     0.01     0.01  0.0  0.     0.01     0.01  0.0 \n')	
			f1.close()
			
			print 'FILE NOT FOUND'
			continue
			
		mcmc_both           = mcmc_both[500:1000] 
		mcmc_both_miss      = mcmc_both_miss[500:1000] 
		mcmc_both_in        = mcmc_both_in[500:1000] 
		mcmc_both_miss_in   = mcmc_both_miss_in[500:1000] 
		mcmc_both_out       = mcmc_both_out[500:1000] 
		mcmc_both_miss_out  = mcmc_both_miss_out[500:1000] 
		
		eb          = np.percentile(mcmc_both, [16, 50, 84])
		eb_miss     = np.percentile(mcmc_both_miss, [16, 50, 84])
		eb_in       = np.percentile(mcmc_both_in, [16, 50, 84])
		eb_miss_in  = np.percentile(mcmc_both_miss_in, [16, 50, 84])
		eb_out      = np.percentile(mcmc_both_out, [16, 50, 84])
		eb_miss_out = np.percentile(mcmc_both_miss_out, [16, 50, 84])
		
		e = [eb[1],eb_miss[1],eb_in[1],
		     eb_miss_in[1],eb_out[1],eb_miss_out[1]]
		
		f, ax = plt.subplots(figsize=(3,5))
		ax.hist(mcmc_both,histtype='step',alpha=0.7,color='C0',label='centred')
		ax.hist(mcmc_both_miss,histtype='step',alpha=0.7,color='C1',label='misscentred')
		ax.hist(mcmc_both_in,histtype='step',alpha=0.7,color='C2',label='centred_in')
		ax.hist(mcmc_both_miss_in,histtype='step',alpha=0.7,color='C3',label='misscentred_in')
		ax.hist(mcmc_both_out,histtype='step',alpha=0.7,color='C4',label='centred_out')
		ax.hist(mcmc_both_miss_out,histtype='step',alpha=0.7,color='C5',label='misscentred_out')
		ax.set_xlabel('e')
		f.subplots_adjust(hspace=0,wspace=0)
		plt.savefig(folder+'plots_all_mcmc/hist_'+plot_out)
		
		
		profile = np.loadtxt(folder+file_profile).T
		
		f = open(folder+file_profile,'r')
		lines = f.readlines()
		j = lines[1].find('=')+1
		Mguess = (float(lines[1][j:-2])*1.e14)
		j = lines[2].find('=')+1
		zmean = float(lines[2][j:-2])
	
		######################
		#COMPUTE CHI_SQUARE
	
		chi = []
		miss = [False,True,False,True,False,True]
		r  = np.logspace(np.log10(min(profile[0])),
						np.log10(max(profile[0])),20)

		ftype = ['C2','C4','C2--','C4--','C2-.','C4-.']

		f, ax = plt.subplots(1, 2, figsize=(8,5))
		ax[0].scatter(profile[0],profile[1],facecolor='none',edgecolors='0.4')
		ax[0].errorbar(profile[0],profile[1],yerr=profile[2],fmt = 'none',ecolor='0.4')
		ax[1].scatter(profile[0],profile[3],facecolor='none',edgecolors='0.4')
		ax[1].errorbar(profile[0],profile[3],yerr=profile[4],fmt = 'none',ecolor='0.4')


		
		for j in range(6):
	
			multipoles = multipole_shear_parallel(profile[0],M200=M200,misscentred = miss[j],
										ellip=e[j],z=zmean,components = ['tcos','xsin'],
										verbose=False,ncores=ncores)
			
			modelb_t = model_Gamma(multipoles,'tcos', misscentred = miss[j], pcc = pcc)
			modelb_x = model_Gamma(multipoles,'xsin', misscentred = miss[j], pcc = pcc)
		
			chi = np.append(chi,chi_red(np.append(modelb_t,modelb_x),
					np.append(profile[1],profile[3]),
					np.append(profile[2],profile[4]),1))

			multipoles = multipole_shear_parallel(r,M200=M200,misscentred = miss[j],
										ellip=e[j],z=zmean,components = ['tcos','xsin'],
										verbose=False,ncores=ncores)
			
			modelb_t = model_Gamma(multipoles,'tcos', misscentred = miss[j], pcc = pcc)
			modelb_x = model_Gamma(multipoles,'xsin', misscentred = miss[j], pcc = pcc)

			ax[0].plot(r,modelb_t,ftype[j])
			ax[1].plot(r,modelb_x,ftype[j])
		######################		
		
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
		
		ax[1].set_xlabel('R [mpc]')
		ax[1].set_xscale('log')
		ax[1].set_ylim(-50,50)
		ax[1].set_xlim(0.1,5.2)
		ax[1].xaxis.set_ticks([0.1,1,5])
		ax[1].set_xticklabels([0.1,1,5])
		ax[1].yaxis.set_ticks([-30,-15,0,15,30])
		ax[1].set_yticklabels([-30,-15,0,15,30])
		plt.legend()
		f.subplots_adjust(hspace=0,wspace=0)
		plt.savefig(folder+'plots_all_mcmc/'+plot_out)
		
		#/////// save file ///////
		
		f1=open(folder+out_file,'a')
		f1.write(file_profile+' ') 
		f1.write(str('%.2f' % (M200/1.e14))+'   '+str('%.2f' % (zmean))+'   0   5    ')
		f1.write(str('%.2f' % (eb[1]))+'   '+str('%.2f' % (np.diff(eb)[0]))+'   '+str('%.2f' % (np.diff(eb)[1]))+'   ')
		f1.write(str('%.2f' % chi[0])+'   ')
		f1.write(str('%.2f' % (eb_miss[1]))+'   '+str('%.2f' % (np.diff(eb_miss)[0]))+'   '+str('%.2f' % (np.diff(eb_miss)[1]))+'   ')
		f1.write(str('%.2f' % chi[1])+'   \n')

		f1.write(file_profile+' ') 
		f1.write(str('%.2f' % (M200/1.e14))+'   '+str('%.2f' % (zmean))+'   0   1    ')
		f1.write(str('%.2f' % (eb_in[1]))+'   '+str('%.2f' % (np.diff(eb_in)[0]))+'   '+str('%.2f' % (np.diff(eb_in)[1]))+'   ')
		f1.write(str('%.2f' % chi[2])+'   ')
		f1.write(str('%.2f' % (eb_miss_in[1]))+'   '+str('%.2f' % (np.diff(eb_miss_in)[0]))+'   '+str('%.2f' % (np.diff(eb_miss_in)[1]))+'   ')
		f1.write(str('%.2f' % chi[3])+'   \n')

		f1.write(file_profile+' ') 
		f1.write(str('%.2f' % (M200/1.e14))+'   '+str('%.2f' % (zmean))+'   1   5    ')
		f1.write(str('%.2f' % (eb_out[1]))+'   '+str('%.2f' % (np.diff(eb_out)[0]))+'   '+str('%.2f' % (np.diff(eb_out)[1]))+'   ')
		f1.write(str('%.2f' % chi[4])+'   ')
		f1.write(str('%.2f' % (eb_miss_out[1]))+'   '+str('%.2f' % (np.diff(eb_miss_out)[0]))+'   '+str('%.2f' % (np.diff(eb_miss_out)[1]))+'   ')
		f1.write(str('%.2f' % chi[5])+'   \n')

		f1.close()
		print 'tardo... '
		print (time.time() -t1)/60.


def final_plots(folder,file_name,ang,miss,
				 ncores,ax1,ax2,prox):
	
	print ang
	file_profile = file_name[:-4]+'_'+ang+'.cat'
	print file_profile

	if miss:
		file_mcmc_in  = folder+'quadrupole_both_miss_'+file_name[:-4]+'_'+ang+'_0_700.out'
		file_mcmc_out = folder+'quadrupole_both_miss_'+file_name[:-4]+'_'+ang+'_700_5000.out'
		file_mcmc_tot = folder+'quadrupole_both_miss_'+file_name[:-4]+'_'+ang+'_0_5000.out'
	else:
		file_mcmc_in  = folder+'quadrupole_both_'+file_name[:-4]+'_'+ang+'_0_700.out'
		file_mcmc_out = folder+'quadrupole_both_'+file_name[:-4]+'_'+ang+'_700_5000.out'
		file_mcmc_tot = folder+'quadrupole_both_'+file_name[:-4]+'_'+ang+'_0_5000.out'

	
	f = open(folder+file_name,'r')
	lines = f.readlines()
	j = lines[2].find('=')+1
	zmean = float(lines[2][j:-2])
	pcc = float((lines[-1][1:-2]))
	M200 = float((lines[-2][1:-2]))*1.e14
	
	mcmc_in  = (np.loadtxt(file_mcmc_in))[500:]
	mcmc_out = (np.loadtxt(file_mcmc_out))[500:]
	mcmc_tot = (np.loadtxt(file_mcmc_tot))[500:]
	
	
	e_in   = np.percentile(mcmc_in,  [16, 50, 84])
	e_out  = np.percentile(mcmc_out, [16, 50, 84])
	e_tot  = np.percentile(mcmc_tot, [16, 50, 84])
		
	profile = np.loadtxt(folder+file_profile).T
	
	f = open(folder+file_profile,'r')
	lines = f.readlines()
	j = lines[1].find('=')+1
	Mguess = (float(lines[1][j:-2])*1.e14)
	j = lines[2].find('=')+1
	zmean = float(lines[2][j:-2])

	
	r  = np.logspace(np.log10(0.01),np.log10(5.5),20)
	
	multipoles = multipole_shear_parallel(r,M200=M200,misscentred = miss,
								ellip=e_in[1],z=zmean,components = ['tcos','xsin'],
								verbose=False,ncores=ncores)
	
	modelin_t = model_Gamma(multipoles,'tcos', misscentred = miss, pcc = pcc)
	modelin_x = model_Gamma(multipoles,'xsin', misscentred = miss, pcc = pcc)

	multipoles = multipole_shear_parallel(r,M200=M200,misscentred = miss,
								ellip=e_out[1],z=zmean,components = ['tcos','xsin'],
								verbose=False,ncores=ncores)
	
	modelout_t = model_Gamma(multipoles,'tcos', misscentred = miss, pcc = pcc)
	modelout_x = model_Gamma(multipoles,'xsin', misscentred = miss, pcc = pcc)

	multipoles = multipole_shear_parallel(r,M200=M200,misscentred = miss,
								ellip=e_tot[1],z=zmean,components = ['tcos','xsin'],
								verbose=False,ncores=ncores)
	
	modeltot_t = model_Gamma(multipoles,'tcos', misscentred = miss, pcc = pcc)
	modeltot_x = model_Gamma(multipoles,'xsin', misscentred = miss, pcc = pcc)
	
	ax1.scatter(profile[0],profile[1],facecolor='none',edgecolors='0.4')
	ax1.plot(r,modeltot_t,'C3',label = r'$r^{0.1}_{5.0}$')
	ax1.fill_between(r, modeltot_t - modeltot_t*np.diff(e_tot)[0], modeltot_t + modeltot_t*np.diff(e_tot)[1],facecolor = 'C3', alpha = 0.3)
	ax1.plot(r,modelin_t,'C4',label = r'$r^{0.1}_{0.7}$')
	ax1.fill_between(r, modelin_t - modelin_t*np.diff(e_in)[0], modelin_t + modelin_t*np.diff(e_in)[1],facecolor = 'C4', alpha = 0.3)
	ax1.plot(r,modelout_t,'C5',label = r'$r^{0.7}_{5.0}$')
	ax1.fill_between(r, modelout_t - modelout_t*np.diff(e_out)[0], modelout_t + modelout_t*np.diff(e_out)[1],facecolor = 'C5', alpha = 0.3)
	ax1.errorbar(profile[0],profile[1],yerr=profile[2],fmt = 'none',ecolor='0.4')
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	# ax1.set_xlabel(r'$r [h^{-1}_{70}\,Mpc]$')
	ax1.set_ylim(0.1,200)
	ax1.set_xlim(0.1,5.2)
	ax1.xaxis.set_ticks([0.2,1,4])
	ax1.set_xticklabels([0.2,1,4])
	ax1.text(2,60,prox)
	ax1.yaxis.set_ticks([1,10,100])
	ax1.set_yticklabels([1,10,100])

	
	ax2.scatter(profile[0],profile[3],facecolor='none',edgecolors='0.4')
	ax2.plot([0,5],[0,0],'C7--')
	ax2.plot(r,modeltot_x,'C3',label = r'$r^{0.1}_{5.0}$')
	ax2.fill_between(r, modeltot_x - modeltot_x*np.diff(e_tot)[0], modeltot_x + modeltot_x*np.diff(e_tot)[1],facecolor = 'C3', alpha = 0.3)
	ax2.plot(r,modelin_x,'C4',label = r'$r^{0.1}_{0.7}$')
	ax2.fill_between(r, modelin_x - modelin_x*np.diff(e_in)[0], modelin_x + modelin_x*np.diff(e_in)[1],facecolor = 'C4', alpha = 0.3)
	ax2.plot(r,modelout_x,'C5',label = r'$r^{0.7}_{5.0}$')
	ax2.fill_between(r, modelout_x - modelout_x*np.diff(e_out)[0], modelout_x + modelout_x*np.diff(e_out)[1],facecolor = 'C5', alpha = 0.3)
	ax2.errorbar(profile[0],profile[3],yerr=profile[4],fmt = 'none',ecolor='0.4')
	ax2.text(2,20,prox)
	ax2.set_xlabel(r'$r [h^{-1}_{70}\,Mpc]$')
	ax2.set_xscale('log')
	ax2.set_ylim(-70,50)
	ax2.set_xlim(0.1,5.2)
	ax2.xaxis.set_ticks([0.2,1,4])
	ax2.set_xticklabels([0.2,1,4])
	ax2.yaxis.set_ticks([-60,-30,0,30])
	ax2.set_yticklabels([-60,-30,0,30])

	# f.subplots_adjust(hspace=0,wspace=0)
	

def finalplot2():
	folder = u'/home/eli/Documentos/Astronomia/posdoc/halo-elongation/redMapper/profiles/'
	file_name = 'profile_total.cat'
	quadrupole = np.loadtxt(folder+'profile_total_t.cat').T

	f = open(folder+file_name,'r')
	lines = f.readlines()
	j = lines[2].find('=')+1
	zmean = float(lines[2][j:-2])
	pcc = float((lines[-1][1:-2]))
	M200 = float((lines[-2][1:-2]))*1.e14

	e = 0.17
	r  = np.logspace(np.log10(0.01),np.log10(5.5),20)
	
	multipoles = multipole_shear_parallel(r,M200=M200,misscentred = True, ellip=0,z=zmean,components = ['t'],ncores=2)
	mono = model_Gamma(multipoles,'t', misscentred = True, pcc = pcc)
	mono_cen = multipoles['Gt0']
	mono_miss = multipoles['Gt_off']
	
	multipoles = multipole_shear_parallel(r,M200=M200,misscentred = False, ellip=e,z=zmean,ncores=2)

	f, ax = plt.subplots(3,1, figsize=(3.5,7), sharex=True)
	f.subplots_adjust(hspace=0,wspace=0)
	matplotlib.rcParams.update({'font.size': 12})
	
	
	ax[0].scatter(profile[0],profile[1],facecolor='none',edgecolors='0.4')
	ax[0].errorbar(profile[0],profile[1],yerr=profile[2],fmt = 'none',ecolor='0.4')
	ax[0].plot(r,mono,'C3')
	ax[0].plot(r,pcc*mono_cen,'C4')
	ax[0].plot(r,(1-pcc)*mono_miss,'C4--')
	ax[0].set_xscale('log')
	ax[0].set_yscale('log')
	ax[0].set_ylabel(r'$\Delta \Sigma [h_{70}M_\odot/pc]$')
	
	ax[1].scatter(quadrupole[0],quadrupole[1],facecolor='none',edgecolors='0.4')
	ax[1].errorbar(quadrupole[0],quadrupole[1],yerr=quadrupole[2],fmt = 'none',ecolor='0.4')
	ax[1].plot(r,multipoles['Gt2'],'C3')
	ax[1].set_xscale('log')
	ax[1].set_yscale('log')
	ax[1].set_ylabel(r'$\Gamma_{\rm{t} \cos{2\theta}} [h_{70}M_\odot/pc]$')

	plt.rc('font', family='serif', size='12.0')
	ax[2].scatter(quadrupole[0],quadrupole[3],facecolor='none',edgecolors='0.4')
	ax[2].errorbar(quadrupole[0],quadrupole[3],yerr=quadrupole[4],fmt = 'none',ecolor='0.4')
	ax[2].plot(r,multipoles['Gx2'],'C3')
	ax[2].plot([0,6],[0,0],'C7--')
	ax[2].set_ylim(-70,30)
	ax[2].yaxis.set_ticks([-60,-30,0])
	ax[2].set_yticklabels([-60,-30,0])
	ax[2].set_ylabel(r'$\Gamma_{\rm{\times} \sin{2\theta}} [h_{70}M_\odot/pc]$') 
	ax[2].set_xlabel(r'$r [h^{-1}_{70}\,Mpc]$')
	
	ax[0].xaxis.set_ticks([0.2,1,4])
	ax[0].set_xticklabels([0.2,1,4])
	ax[0].set_ylim(1,300)
	ax[0].set_xlim(0.1,5.2)
	
	f.savefig('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/redMapper/profiles.pdf',bbox_inches='tight')


