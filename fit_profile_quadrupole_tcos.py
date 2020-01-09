import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
sys.path.append('/home/elizabeth/multipole_density_profile')
sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
import numpy as np
from pylab import *
from multipoles_shear import *
import emcee
import time
from multiprocessing import Pool
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-folder', action='store', dest='folder',default='./')
parser.add_argument('-file', action='store', dest='file_name', default='profile.cat')
parser.add_argument('-ncores', action='store', dest='ncores', default=4)
args = parser.parse_args()

folder    = args.folder
file_name = args.file_name
ncores    = args.ncores
ncores    = int(ncores)

print 'fitting quadrupole tcos'
print folder
print file_name


f = open(folder+file_name,'r')
lines = f.readlines()
j = lines[1].find('=')+1
Mguess = (float(lines[1][j:-2])*1.e14)*1.3
j = lines[2].find('=')+1
zmean = float(lines[2][j:-2])
    

def log_likelihood(data_model, r, Gamma, e_Gamma):
    log_M200,ellip = data_model
    M200 = 10**log_M200
    multipoles = multipole_shear(r,M200=M200,misscentred = False,
                                ellip=ellip,z=zmean,components = ['tcos'],
                                verbose=False)
    model = model_Gamma(multipoles,'tcos', misscentred = False)
    sigma2 = e_Gamma**2
    return -0.5 * np.sum((Gamma - model)**2 / sigma2 + np.log(2.*np.pi*sigma2))
    

def log_probability(data_model, r, Gamma, e_Gamma):
    log_M200, ellip = data_model
    if np.log10(Mguess*0.5) < log_M200 < np.log10(Mguess*1.5) and 0. < ellip < 0.5:
        return log_likelihood(data_model, r, Gamma, e_Gamma)
    return -np.inf

# initializing

pos = np.array([np.random.normal(np.log10(Mguess),0.2,30),
                np.random.uniform(0.,0.5,30)]).T


nwalkers, ndim = pos.shape

#-------------------
# running emcee

pool = Pool(processes=(ncores))

profile = np.loadtxt(folder+file_name).T

t1 = time.time()
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, 
                                args=(profile[0],profile[1],profile[2]),
                                pool = pool)
sampler.run_mcmc(pos, 1000, progress=True)
print (time.time()-t1)/60.
pool.terminate()
#-------------------

#flat_samples = sampler.get_chain(discard=100, flat=True)

#p1 = np.percentile(flat_samples[:, 0], [16, 50, 84])

#print p1[1],np.diff(p1)

# saving mcmc out

mcmc_out = sampler.get_chain(flat=True)

f1=open(folder+'quadrupole_tcos_'+file_name,'w')
f1.write('# log(M200)  ellip \n')
np.savetxt(f1,mcmc_out,fmt = ['%12.6f']*2)
f1.close()
