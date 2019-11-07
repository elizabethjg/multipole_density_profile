import numpy as np
from pylab import *
from astropy.cosmology import LambdaCDM
from scipy.misc import derivative
from scipy import integrate
from multipoles_shear import *
from multiprocessing import Pool
from multiprocessing import Process

#parameters

cvel = 299792458;   # Speed of light (m.s-1)
G    = 6.670e-11;   # Gravitational constant (m3.kg-1.s-2)
pc   = 3.085678e16; # 1 pc (m)
Msun = 1.989e30 # Solar mass (kg)

def pru(r,a=2):
     return {'r':r, '2a':a*2.}

def pru2(into):
     return pru(*into)

ncores = 2

def run_parallel(ncores=2):
     r = np.array([0.1,0.2,0.3,0.4])
     slicer = int(round(len(r)/ncores, 0))
     slices = ((np.arange(ncores-1)+1)*slicer).astype(int)
     r_splited = np.split(r,slices)
     a = np.ones(len(r_splited))
     
     pool = Pool(processes=(ncores))
     salida=np.array(pool.map(pru2, np.array([r_splited,a]).T))
     pool.terminate()
     return salida

# '''
############  MAKING A GRID

logM = 14.
M    = 10**logM
q    = 0.6

a  = np.logspace(np.log10(0.05),np.log10(1.5),10)
a  = np.append(a,-1.*a)

x,y = np.meshgrid(a,a)

x = x.flatten()
y = y.flatten()

r = np.sqrt(x**2 + y**2)

mask = (r>0.01)

x,y = x[mask],y[mask]
r   = r[mask]
theta  = np.arctan2(y,x)

j   = argsort(r)
r   = r[j]
theta  = theta[j]
x,y = x[j],y[j]




'''
