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
parser.add_argument('-misscentred', action='store', dest='miss', default=0)
parser.add_argument('-component', action='store', dest='component', default='tcos')
parser.add_argument('-RIN', action='store', dest='RIN', default=0)
parser.add_argument('-ROUT', action='store', dest='ROUT', default=5)
args = parser.parse_args()

folder    = args.folder
file_name = args.file_name
angle     = args.angle

if 'True' in args.miss:
	miss      = True
elif 'False' in args.miss:
	miss      = False
	
component = args.component
ncores    = args.ncores
ncores    = int(ncores)
rin       = args.RIN
rout      = args.ROUT

print 'fitting quadrupole'
print folder
print file_name
print angle
print 'miss = ',miss
print component
print 'ncores = ',ncores
print 'RIN ',rin
print 'ROUT ',rout


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
    if 'both' in component:
	r = np.split(r,2)[0]
	multipoles = multipole_shear_parallel(r,M200=M200,misscentred = miss,
				    ellip=ellip,z=zmean,components = ['tcos','xsin'],
				    verbose=False,ncores=ncores)
	model_t = model_Gamma(multipoles,'tcos', misscentred = miss, pcc = pcc)
	model_x = model_Gamma(multipoles,'xsin', misscentred = miss, pcc = pcc)
	model   = np.append(model_t,model_x)
    else:
	multipoles = multipole_shear_parallel(r,M200=M200,misscentred = miss,
				    ellip=ellip,z=zmean,components = [component],
				    verbose=False,ncores=ncores)
	model = model_Gamma(multipoles,component, misscentred = miss, pcc = pcc)
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
print rin,rout
maskr   = (profile[0]>rin)*(profile[0]<rout)
print profile[0]
print maskr.sum()
print (profile[0]>rin)
print (profile[0]<rout)
print type(rin)
print type(rout)
profile = profile[:,maskr]

t1 = time.time()

if component == 'tcos':
	
	sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, 
                                args=(profile[0],profile[1],profile[2]))
elif component == 'xsin':                                
	sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, 
                                args=(profile[0],profile[3],profile[4]))
elif component == 'both':                                
	sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, 
                                args=(np.append(profile[0],profile[0]),
				np.append(profile[1],profile[3]),
				np.append(profile[2],profile[4])))


                                
sampler.run_mcmc(pos, 250, progress=True)
print (time.time()-t1)/60.

#-------------------
# saving mcmc out

mcmc_out = sampler.get_chain(flat=True)

if miss:
	f1=open(folder+'quadrupole_'+component+'_miss_'+file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'.out','w')
else:
	f1=open(folder+'quadrupole_'+component+'_'+file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'.out','w')
f1.write('# ellip \n')
np.savetxt(f1,mcmc_out,fmt = ['%12.6f'])
f1.close()
