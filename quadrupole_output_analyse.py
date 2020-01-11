import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
sys.path.append('/home/elizabeth/multipole_density_profile')
sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
import numpy as np
from pylab import *
from multipoles_shear import *
import time
import corner
import os

folder         = u'/home/eli/Documentos/PostDoc/halo-elongation/redMapper/profiles_newbins/'
file_profile   = 'profile_bin3_tpwl.cat'
out_file       = 'out_mcmc_table'

file_mcmc_xsin = 'quadrupole_xsin_'+file_profile
file_mcmc_tcos = 'quadrupole_tcos_'+file_profile


os.system('mkdir '+folder+'plots_mcmc/')

mcmc_xsin = (np.loadtxt(folder+file_mcmc_xsin)[3000:]).T
mcmc_tcos = (np.loadtxt(folder+file_mcmc_tcos)[3000:]).T
mcmc_combined = np.concatenate((mcmc_tcos,mcmc_tcos),axis=1)


labels = ['M200','e']

#/////// xsin ///////

m200_x = np.percentile(mcmc_xsin[0], [16, 50, 84])
e_x    = np.percentile(mcmc_xsin[1], [16, 50, 84])

m200x   = 10**(m200_x[1])/1.e14
e_m200x = (10**(m200_x[1])*np.log(10.)*np.diff(m200_x))/1.e14

fig = corner.corner(mcmc_xsin.T, labels=labels)
plt.savefig(folder+'plots_mcmc/'+file_mcmc_xsin[:-3]+'png')

#/////// tcos ///////

m200_t = np.percentile(mcmc_tcos[0], [16, 50, 84])
e_t    = np.percentile(mcmc_tcos[1], [16, 50, 84])

m200t   = 10**(m200_t[1])/1.e14
e_m200t = (10**(m200_t[1])*np.log(10.)*np.diff(m200_t))/1.e14

fig = corner.corner(mcmc_tcos.T)
plt.savefig(folder+'plots_mcmc/'+file_mcmc_xsin[:-3]+'png')

#/////// both ///////

m200_total   = np.percentile(mcmc_combined[0], [16, 50, 84])
e      = np.percentile(mcmc_combined[1], [16, 50, 84])

m200   = 10**(m200_total[1])/1.e14
e_m200 = (10**(m200_total[1])*np.log(10.)*np.diff(m200_total))/1.e14

fig = corner.corner(mcmc_combined.T)
plt.savefig(folder+'plots_mcmc/'+file_mcmc_xsin[:-3]+'png')

#/////// save file ///////

f1=open(folder+out_file,'a')
f1.write(file_profile+'  '+'  tcos  ') 
f1.write(str('%.2f' % (m200t))+'   '+str('%.2f' % (e_m200t[0]))+'   '+str('%.2f' % (e_m200t[1]))+'   ')
f1.write(str('%.2f' % (e_t[1]))+'   '+str('%.2f' % (np.diff(e_t)[0]))+'   '+str('%.2f' % (np.diff(e_t)[1]))+'   ')
f1.write('  xsin  ') 
f1.write(str('%.2f' % (m200x))+'   '+str('%.2f' % (e_m200x[0]))+'   '+str('%.2f' % (e_m200x[1]))+'   ')
f1.write(str('%.2f' % (e_x[1]))+'   '+str('%.2f' % (np.diff(e_x)[0]))+'   '+str('%.2f' % (np.diff(e_x)[1]))+'   ')
f1.write('  both  ') 
f1.write(str('%.2f' % (m200))+'   '+str('%.2f' % (e_m200[0]))+'   '+str('%.2f' % (e_m200[1]))+'   ')
f1.write(str('%.2f' % (e[1]))+'   '+str('%.2f' % (np.diff(e)[0]))+'   '+str('%.2f' % (np.diff(e)[1]))+'   \n')
f1.close()


profile = np.loadtxt(folder+file_profile).T

f = open(folder+file_profile,'r')
lines = f.readlines()
j = lines[1].find('=')+1
Mguess = (float(lines[1][j:-2])*1.e14)
j = lines[2].find('=')+1
zmean = float(lines[2][j:-2])

r  = np.logspace(np.log10(min(profile[0])),
				np.log10(max(profile[0])),10)

multipoles = multipole_shear_parallel(r,M200=m200t*1.e14,components = ['tcos','xsin'],
									misscentred = False,
									ellip=e_t[1],z=zmean,
									verbose=True,ncores=4)

modelt_t = model_Gamma(multipoles,'tcos', misscentred = False)
modelt_x = model_Gamma(multipoles,'xsin', misscentred = False)

multipoles = multipole_shear_parallel(r,M200=m200x*1.e14,components = ['tcos','xsin'],
									misscentred = False,
									ellip=e_x[1],z=zmean,
									verbose=True,ncores=4)

modelx_t = model_Gamma(multipoles,'tcos', misscentred = False)
modelx_x = model_Gamma(multipoles,'xsin', misscentred = False)


multipoles = multipole_shear_parallel(r,M200=m200*1.e14,components = ['tcos','xsin'],
									misscentred = False,
									ellip=e[1],z=zmean,
									verbose=True,ncores=4)

model_t = model_Gamma(multipoles,'tcos', misscentred = False)
model_x = model_Gamma(multipoles,'xsin', misscentred = False)

f, ax = plt.subplots(figsize=(6.5,5))
ax.plot(profile[0],profile[1],'C0o')
ax.plot(r,model_t,'C1',label = 'combined')
ax.plot(r,modelt_t,'C2',label = 'tcos')
ax.plot(r,modelx_t,'C3',label = 'xsin')
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
plt.legend()

f, ax = plt.subplots(figsize=(6.5,5))
ax.plot(profile[0],profile[3],'C0o')
ax.plot(r,model_x,'C1',label = 'combined')
ax.plot(r,modelt_x,'C2',label = 'tcos')
ax.plot(r,modelx_x,'C3',label = 'xsin')
ax.errorbar(profile[0],profile[3],yerr=profile[4],fmt = 'none',ecolor='C0')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('R [mpc]')
ax.set_ylim(1,300)
ax.set_xlim(0.1,5.2)
ax.xaxis.set_ticks([0.1,1,5])
ax.set_xticklabels([0.1,1,5])
ax.yaxis.set_ticks([1,10,100])
ax.set_yticklabels([1,10,100])
plt.legend()
