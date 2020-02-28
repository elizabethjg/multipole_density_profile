import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
sys.path.append('/home/elizabeth/multipole_density_profile')
sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
import numpy as np
import matplotlib
matplotlib.use('Agg')
from pylab import *
from multipoles_shear import *
import time
import os
from profiles_fit import chi_red

def plot_mcmc_quadrupole_out(folder,file_name,angle,miss,
							 rin,rout,out_file,ncores):
	
	file_profile = file_name[:-4]+'_'+angle+'.cat'
	
	print file_profile

	if miss:
		file_mcmc_xsin = folder+'quadrupole_xsin_miss_'+file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		file_mcmc_tcos = folder+'quadrupole_tcos_miss_'+file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		file_mcmc_both = folder+'quadrupole_both_miss_'+file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		plot_out       = file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'_miss.png'
	else:
		file_mcmc_xsin = folder+'quadrupole_xsin_'+file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		file_mcmc_tcos = folder+'quadrupole_tcos_'+file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		file_mcmc_both = folder+'quadrupole_both_'+file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		plot_out       = file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'.png'
	
	os.system('mkdir '+folder+'plots_mcmc/')
	
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
		return None
	
	
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
	
	f, ax = plt.subplots(1, 3, figsize=(10,5))
	
	ax[0].hist(mcmc_tcos,histtype='step')
	ax[0].set_xlabel('e_t')
	ax[1].hist(mcmc_xsin,histtype='step')
	ax[1].set_xlabel('e_x')
	ax[2].hist(mcmc_both,histtype='step')
	ax[2].set_xlabel('e_b')
	f.subplots_adjust(hspace=0,wspace=0)
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
	ax[0].plot(r,modelx_t,'C3',label = 'xsin')
	ax[0].plot(r,modelb_t,'C4',label = 'both')
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
	ax[1].plot(r,modelx_x,'C3',label = 'xsin')
	ax[1].plot(r,modelb_x,'C4',label = 'both')
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
	
	#/////// save file ///////
	
	f1=open(folder+out_file,'a')
	f1.write(file_profile+' ') 
	f1.write(str('%.2f' % (M200/1.e14))+'   '+str('%.2f' % (zmean))+'   '+str(int(rin))+'   '+str(int(rout))+'    ')
	f1.write(str('%.2f' % (e_t[1]))+'   '+str('%.2f' % (np.diff(e_t)[0]))+'   '+str('%.2f' % (np.diff(e_t)[1]))+'   '+str('%.2f' % (chi_t))+'   ')
	f1.write(str('%.2f' % (e_x[1]))+'   '+str('%.2f' % (np.diff(e_x)[0]))+'   '+str('%.2f' % (np.diff(e_x)[1]))+'   '+str('%.2f' % (chi_x))+'   ')
	f1.write(str('%.2f' % (e_b[1]))+'   '+str('%.2f' % (np.diff(e_b)[0]))+'   '+str('%.2f' % (np.diff(e_b)[1]))+'   '+str('%.2f' % (chi_b))+' \n')
	f1.close()

def plot_mcmc_quadrupole_out_onlyboth(folder,file_name,angle,miss,
							 rin,rout,out_file,ncores):
	
	file_profile = file_name[:-4]+'_'+angle+'.cat'
	
	print file_profile

	if miss:
		file_mcmc_xsin = folder+'quadrupole_xsin_miss_'+file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		file_mcmc_tcos = folder+'quadrupole_tcos_miss_'+file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		file_mcmc_both = folder+'quadrupole_both_miss_'+file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		plot_out       = file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'_miss.png'
	else:
		file_mcmc_xsin = folder+'quadrupole_xsin_'+file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		file_mcmc_tcos = folder+'quadrupole_tcos_'+file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		file_mcmc_both = folder+'quadrupole_both_'+file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'.out'
		plot_out       = file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'.png'
	
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
	
		
	
	angles = ['t','twl','twd','tp','tpwl',
              'tpwd','control']
	
	for angle in angles:
		
		t1 = time.time()
		
		file_profile = file_name[:-4]+'_'+angle+'.cat'
		
		print file_profile

		file_mcmc_both_miss = folder+'quadrupole_both_miss_'+file_name[:-4]+'_'+angle+'_0_5.out'
		file_mcmc_both = folder+'quadrupole_both_'+file_name[:-4]+'_'+angle+'_0_5.out'
		file_mcmc_both_miss_in = folder+'quadrupole_both_miss_'+file_name[:-4]+'_'+angle+'_0_1.out'
		file_mcmc_both_in = folder+'quadrupole_both_'+file_name[:-4]+'_'+angle+'_0_1.out'
		file_mcmc_both_miss_out = folder+'quadrupole_both_miss_'+file_name[:-4]+'_'+angle+'_1_5.out'
		file_mcmc_both_out = folder+'quadrupole_both_'+file_name[:-4]+'_'+angle+'_1_5.out'
		plot_out       = file_name[:-4]+'_'+angle+'.png'
		
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
