import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
sys.path.append('/home/elizabeth/multipole_density_profile')
sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
import numpy as np
from pylab import *
from multipoles_shear import *
import time
from scipy.optimize import curve_fit
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-folder', action='store', dest='folder',default='./')
parser.add_argument('-file', action='store', dest='file_name', default='profile.cat')
parser.add_argument('-ang', action='store', dest='angle', default='twl')
parser.add_argument('-ncores', action='store', dest='ncores', default=4)
parser.add_argument('-misscentred', action='store', dest='miss', default=0)
parser.add_argument('-component', action='store', dest='component', default='both')
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
nit       = int(args.nit)
ncores    = args.ncores
ncores    = int(ncores)
rin       = float(args.RIN)
rout      = float(args.ROUT)

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


profile = np.loadtxt(folder+file_name[:-4]+'_'+angle+'.cat').T
maskr   = (profile[0]>rin)*(profile[0]<rout)
profile = profile[:,maskr]


infit = output = {'r':profile[0],'M200':M200,
                  'z':zmean,'h':0.7,'misscentred':miss,
		  's_off':0.4,'components':[component],
		  'verbose':True,'ncores':ncores}

#-------------------
# running emcee



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


                                
sampler.run_mcmc(pos, nit, progress=True)
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
