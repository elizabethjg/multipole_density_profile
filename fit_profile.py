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

def log_likelihood(data_model, r, Gamma, e_Gamma):
    log_M200, ellip = data_model
    # ellip = (1.- q)/(1. + q)
    M200 = 10**log_M200
    # model = multipole_shear_parallel(r,M200=M200,ellip=ellip,z=0.2,zs=0.6,verbose=False)['Gt2']
    multipoles = multipole_shear(r,M200=M200,ellip=ellip,z=0.2,zs=0.6,verbose=False)
    model = model_Gamma(multipoles,'tcos')
    sigma2 = e_Gamma**2
    return -0.5 * np.sum((Gamma - model)**2 / sigma2 + np.log(2.*np.pi*sigma2))
    

def log_probability(data_model, r, Gamma, e_Gamma):
    log_M200, ellip = data_model
    if 13.5 < log_M200 < 14.5 and 0.0 < ellip < 0.4:
        return log_likelihood(data_model, r, Gamma, e_Gamma)
    return -np.inf

pos = np.array([np.random.normal(13.7,1.0,40),np.random.normal(0.3,0.2,40)]).T
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
