# sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
# sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
sys.path.append('/home/elizabeth/multipole_density_profile')
#sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
# sys.path.append('/home/elizabeth/multipole_density_profile')
import numpy as np
from pylab import *
from multipoles_shear import *
import emcee
import time
from multiprocessing import Pool

def fit_profile_monopole_onlymass(file_name,ncores=2,
                                  folder = './'):
    
	# folder = './'
	# file_name = 'profile_original_bin4.cat'
	# ncores    = 30
	
	f = open(folder+file_name,'r')
	lines = f.readlines()
	j = lines[1].find('=')+1
	Mguess = (float(lines[1][j:-2])*1.e14)
	j = lines[2].find('=')+1
	zmean = float(lines[2][j:-2])
		
	
	def log_likelihood(data_model, r, Gamma, e_Gamma):
		log_M200 = data_model
		M200 = 10**log_M200
		multipoles = multipole_shear(r,M200=M200,misscentred = False,
									ellip=0,z=zmean,verbose=False)
		model = model_Gamma(multipoles,'t', misscentred = False)
		sigma2 = e_Gamma**2
		return -0.5 * np.sum((Gamma - model)**2 / sigma2 + np.log(2.*np.pi*sigma2))
		
	
	def log_probability(data_model, r, Gamma, e_Gamma):
		log_M200 = data_model
		if np.log10(Mguess*0.6) < log_M200 < np.log10(Mguess*1.4):
			return log_likelihood(data_model, r, Gamma, e_Gamma)
		return -np.inf

	# initializing
	
	pos = np.array([np.random.normal(np.log10(Mguess),0.1,50)]).T
	nwalkers, ndim = pos.shape
	
	#-------------------
	# running emcee
	
	pool = Pool(processes=(ncores))
	
	profile = np.loadtxt(folder+file_name).T
	
	t1 = time.time()
	sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, 
									args=(profile[0],profile[1],profile[2]),
									pool = pool)
	sampler.run_mcmc(pos, 200, progress=True)
	print (time.time()-t1)/60.
	pool.terminate()
	#-------------------
	
	flat_samples = sampler.get_chain(discard=100, flat=True)
	
	p1 = np.percentile(flat_samples[:, 0], [16, 50, 84])
	
	print p1[1],np.diff(p1)
	
	# saving mcmc out
	
	mcmc_out = sampler.get_chain(flat=True)[:,0]   
	
	f1=open(folder+'monopole_onlymass_mcmc_'+file_name,'w')
	f1.write('# log(M200)  \n')
	np.savetxt(f1,mcmc_out,fmt = ['%12.6f'])
	f1.close()
	

'''	
	# computing fitted profile
	
	M200  = 10**p1[1]
	eM200 = ((10**(p1[1]))*np.log(10.)*np.diff(p1))/1.e14
	print M200,eM200
	
	r  = np.logspace(np.log10(min(profile[0])),
					np.log10(max(profile[0])),10)
	
	multipoles = multipole_shear_parallel(r,M200=M200,
										misscentred = False,
										ellip=0,z=zmean,
										verbose=False,ncores=ncores)
	model = model_Gamma(multipoles,'t', misscentred = False)
	
	
	f1=open(folder+'fitted_'+file_name,'w')
	f1.write('# M200 = '+str('%.2f' % (M200/1.e14))+' \n')
	profile = np.column_stack((profile[0]*1.0e-3,profile[1],profile[5],profile[2],profile[6]))
	np.savetxt(f1,profile,fmt = ['%12.6f']*5)
	f1.close()
	
	
	if plot:
		
			
		f, ax = plt.subplots(figsize=(6.5,5))
		ax.plot(profile[0],profile[1],'C0o')
		#ax.plot(r,model,'C1')
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
'''
def fit_profile_monopole_misscentred(file_name,ncores=2,
                                  folder = './'):
       
    f = open(folder+file_name,'r')
    lines = f.readlines()
    j = lines[1].find('=')+1
    Mguess = (float(lines[1][j:-2])*1.e14)*1.3
    j = lines[2].find('=')+1
    zmean = float(lines[2][j:-2])
        
    
    def log_likelihood(data_model, r, Gamma, e_Gamma):
        log_M200,pcc = data_model
        M200 = 10**log_M200
        multipoles = multipole_shear(r,M200=M200,misscentred = True,
                                    ellip=0,z=zmean,components = ['t'],
                                    verbose=True)
        model = model_Gamma(multipoles,'t', misscentred = True, pcc = pcc)
        sigma2 = e_Gamma**2
        return -0.5 * np.sum((Gamma - model)**2 / sigma2 + np.log(2.*np.pi*sigma2))
        
    
    def log_probability(data_model, r, Gamma, e_Gamma):
        log_M200, pcc = data_model
        if np.log10(Mguess*0.5) < log_M200 < np.log10(Mguess*1.5) and 0.6 < pcc < 0.9:
            return log_likelihood(data_model, r, Gamma, e_Gamma)
        return -np.inf
    
    # initializing
    
    pos = np.array([np.random.normal(np.log10(Mguess),0.2,20),
                    np.random.normal(0.8,0.1,20)]).T
    nwalkers, ndim = pos.shape
    
    #-------------------
    # running emcee
    
    pool = Pool(processes=(ncores))
    
    profile = np.loadtxt(folder+file_name).T
    
    t1 = time.time()
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, 
                                    args=(profile[0],profile[1],profile[2]),
                                    pool = pool)
    sampler.run_mcmc(pos, 100, progress=True)
    print (time.time()-t1)/60.
    pool.terminate()
    #-------------------
    
    flat_samples = sampler.get_chain(discard=100, flat=True)
    
    p1 = np.percentile(flat_samples[:, 0], [16, 50, 84])
    
    print p1[1],np.diff(p1)
    
    # saving mcmc out
    
    mcmc_out = sampler.get_chain(flat=True)[:,0]   
    
    f1=open(folder+'mcmc_m200_'+file_name,'w')
    f1.write('# log(M200)  \n')
    np.savetxt(f1,mcmc_out,fmt = ['%12.6f'])
    f1.close()
    
    mcmc_out = sampler.get_chain(flat=True)[:,1]   
    
    f1=open(folder+'mcmc_pcc_'+file_name,'w')
    f1.write('# pcc  \n')
    np.savetxt(f1,mcmc_out,fmt = ['%12.6f'])
    f1.close()



def fit_profile_multipoles():
    
    file_name = '/home/eli/Documentos/Astronomia/posdoc/halo-elongation/redMapper/profiles_data/profile_bin4.cat'
    f = open(file_name,'r')
    lines = f.readlines()
    j = lines[1].find('=')+1
    Mguess = (float(lines[1][j:-2])*1.e14*1.3)
    j = lines[2].find('=')+1
    zmean = float(lines[2][j:-2])
        
    
    def log_likelihood(data_model, r, Gamma, e_Gamma):
        log_M200, ellip, pcc = data_model
        M200 = 10**log_M200
        multipoles = multipole_shear(r,M200=M200,misscentred = True,ellip=ellip,z=zmean,verbose=False)
        model = model_Gamma(multipoles,'t', misscentred = True, pcc = pcc)
        sigma2 = e_Gamma**2
        return -0.5 * np.sum((Gamma - model)**2 / sigma2 + np.log(2.*np.pi*sigma2))
        
    
    def log_probability(data_model, r, Gamma, e_Gamma):
        log_M200, ellip = data_model
        if np.log10(Mguess*0.6) < log_M200 < np.log10(Mguess*1.4) and 0.0 < ellip < 0.4 and 0.6 < pcc < 0.9:
            return log_likelihood(data_model, r, Gamma, e_Gamma)
        return -np.inf
    
    pos = np.array([np.random.normal(np.log10(Mguess),0.5,50),
                    np.random.normal(0.3,0.2,50),
                    np.random.normal(0.75,0.4,50)]).T
                    
    me = pos[:,1]>1.
    pos[me,1] = pos[me,1] == 1.
    nwalkers, ndim = pos.shape
    
    pool = Pool(processes=(ncores))
    
    t1 = time.time()
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(r,Gt2, e_Gt2))
    sampler.run_mcmc(pos, 500, progress=True)
    print (time.time()-t1)/60.
    
    
    '''
    flat_samples = sampler.get_chain(discard=100, thin=1, flat=True)
    print(flat_samples.shape)
    labels = ['M200','e']
    fig = corner.corner(flat_samples, labels=labels, truths=[14.,0.25])
    
    
    fig, axes = plt.subplots(2, figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)
    
    axes[-1].set_xlabel("step number")
    p1 = np.percentile(flat_samples[:, 0], [16, 50, 84])
    p2 = np.percentile(flat_samples[:, 1], [16, 50, 84])
    
    print p1[1],np.diff(p1)
    print p2[1],np.diff(p2)
    '''
