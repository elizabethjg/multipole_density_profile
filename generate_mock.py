import numpy as np
# from pylab import *
from multipoles_shear import *

#parameters

cvel = 299792458;   # Speed of light (m.s-1)
G    = 6.670e-11;   # Gravitational constant (m3.kg-1.s-2)
pc   = 3.085678e16; # 1 pc (m)
Msun = 1.989e30 # Solar mass (kg)

############  MAKING A GRID

logM = 14.
M    = 10**logM

r  = np.logspace(np.log10(0.05),np.log10(1.5),10)

q = 0.6
     
ellip = (1.- q)/(1. + q)
     
out = multipole_shear_parallel(r,M200=M,z=0.2,
                              ellip=ellip,misscentred=False,
                              ncores=2)

modelo = model_Gamma(out,'tcos',misscentred=False)     
Gt2   = modelo*(np.random.rand(len(r))*np.logspace(np.log10(0.5),np.log10(0.2),10) + np.logspace(np.log10(0.75),np.log10(0.9),10))
e_Gt2 = np.logspace(np.log10(0.5),np.log10(0.2),10)*modelo

plt.plot(r,out['Gt2'],'-')
plt.plot(r,modelo,'-')
plt.errorbar(r,Gt2, yerr=e_Gt2,fmt='none')
plt.plot(r,Gt2, 'o')

