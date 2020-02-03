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
parser.add_argument('-ang', action='store', dest='angle', default='twl')
parser.add_argument('-ncores', action='store', dest='ncores', default=4)
args = parser.parse_args()

folder    = args.folder
file_name = args.file_name
angle     = args.angle
ncores    = args.ncores
ncores    = int(ncores)

print 'fitting quadrupole tcos'
print folder
print file_name


f = open(folder+file_name,'r')
lines = f.readlines()
j = lines[2].find('=')+1
zmean = float(lines[2][j:-2])
pcc = float((lines[-1][1:-2]))
M200 = float((lines[-2][1:-2]))*1.e14
    
print 'M200',M200
print 'pcc',pcc


def log_likelihood(data_model, r, Gamma, e_Gamma):
    ellip = data_model
    multipoles = multipole_shear_parallel(r,M200=M200,misscentred = True,
                                ellip=ellip,z=zmean,components = ['tcos'],
                                verbose=False,ncores=ncores)
    model = model_Gamma(multipoles,'tcos', misscentred = True, pcc = pcc)
    sigma2 = e_Gamma**2
    return -0.5 * np.sum((Gamma - model)**2 / sigma2 + np.log(2.*np.pi*sigma2))
    

def log_probability(data_model, r, Gamma, e_Gamma):
    ellip = data_model
    if 0. < ellip < 0.5:
        return log_likelihood(data_model, r, Gamma, e_Gamma)
    return -np.inf

# initializing

pos = np.array([np.random.uniform(0.,0.5,10)]).T

nwalkers, ndim = pos.shape

#-------------------
# running emcee

profile = np.loadtxt(folder+file_name[:-4]+'_'+angle+'.cat').T

t1 = time.time()
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, 
                                args=(profile[0],profile[1],profile[2]))
sampler.run_mcmc(pos, 200, progress=True)
print (time.time()-t1)/60.

#-------------------
# saving mcmc out

mcmc_out = sampler.get_chain(flat=True)

f1=open(folder+'quadrupole_tcos_'+file_name,'w')
f1.write('# ellip \n')
np.savetxt(f1,mcmc_out,fmt = ['%12.6f'])
f1.close()
