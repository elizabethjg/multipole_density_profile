import numpy as np
from pylab import *
from astropy.cosmology import LambdaCDM
from scipy.misc import derivative
from scipy import integrate
from multipoles_shear import *
from scipy.optimize import minimize
import emcee
import corner
import time
from multiprocessing import Pool
    
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
    sampler.run_mcmc(pos, 500, progress=True);
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
print 'END OF PROGRAM'
