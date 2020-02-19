import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
sys.path.append('/home/elizabeth/multipole_density_profile')
sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
import numpy as np
from pylab import *
from multipoles_shear import *
import os

folder = '/home/eli/Documentos/Astronomia/posdoc/halo-elongation/redMapper/profiles/'

out = np.loadtxt(folder+'out_mcmc_centred',dtype='str').T

samples = out[0]
M200    = out[1].astype(float)
lM200   = np.log10(M200*1.e14)
zmean   = out[2].astype(float)
rin     = out[3].astype(float)
rout    = out[4].astype(float)
et      = out[5].astype(float)
e_et    = np.array([out[6].astype(float),out[7].astype(float)])
ex      = out[8].astype(float)
e_ex    = np.array([out[9].astype(float),out[10].astype(float)])
eb      = out[11].astype(float)
e_eb    = np.array([out[12].astype(float),out[13].astype(float)])

mtotal   = []
mmbin1   = []
mmbin2   = []
mtbin1   = []
mtbin2   = []
mtbin3   = []
mmz1     = []
mmz2     = []
mtz1     = []
mtz2     = []
mtz3     = []
mt       = []
mtwl     = []
mtwd     = []
mtp      = []
mtpwl    = []
mtpwd    = []
# mtpwdl   = []
mcontrol = []

angles = ['t','twl','twd','tp','tpwl',
          'tpwd','control']

ang_ind = np.arange(len(angles))

for j in samples:
	mtotal += ['total' in j]
	mmbin1 += ['median_bin1' in j]
	mmbin2 += ['median_bin2' in j]
	mtbin1 += ['terciles_bin1' in j]
	mtbin2 += ['terciles_bin2' in j]
	mtbin3 += ['terciles_bin3' in j]
	mmz1   += ['medianz1' in j]
	mmz2   += ['medianz2' in j]
	mtz1   += ['terciles_z1' in j]
	mtz2   += ['terciles_z2' in j]
	mtz3   += ['terciles_z3' in j]
	mt       +=  ['t.cat' in j]     
	mtwl     +=  ['twl.cat' in j]
	mtwd     +=  ['twd.cat' in j]    
	mtp      +=  ['tp.cat' in j]
	mtpwl    +=  ['tpwl.cat' in j]
	mtpwd    +=  ['tpwd.cat' in j]
	# mtpwdl   +=  ['tpwdl.cat' in j]
	mcontrol +=  ['control.cat' in j]
	
mint  = (rout == 1)
mext  = (rin == 1)
mrtot = (rin == 0)*(rout == 5)


plt.figure()
plt.plot(ang_ind,et[mtotal*mrtot],c='C4',label = 'lM200_mean ='+str(lM200[mtotal*mrtot][0]))
plt.fill_between(ang_ind, (et-e_et[0,:])[mtotal*mrtot], (et+e_et[1,:])[mtotal*mrtot],facecolor = 'C4', alpha = 0.3)
plt.plot(ang_ind,et[mmbin1*mrtot],c='C5',label = 'lM200_mean ='+str(lM200[mmbin1*mrtot][0]))
plt.fill_between(ang_ind, (et-e_et[0,:])[mmbin1*mrtot], (et+e_et[1,:])[mmbin1*mrtot],facecolor = 'C5', alpha = 0.3)
plt.plot(ang_ind,et[mmbin2*mrtot],c='C6',label = 'lM200_mean ='+str(lM200[mmbin2*mrtot][0]))
plt.fill_between(ang_ind, (et-e_et[0,:])[mmbin2*mrtot], (et+e_et[1,:])[mmbin2*mrtot],facecolor = 'C6', alpha = 0.3)
plt.xticks(ang_ind,angles)
plt.legend()
plt.savefig(folder+'plot_results/Mrtot_t.png')

plt.figure()
plt.plot(ang_ind,et[mtotal*mrtot],c='C4',label = 'lM200_mean ='+str(lM200[mtotal*mrtot][0]))
plt.fill_between(ang_ind, (et-e_et[0,:])[mtotal*mrtot], (et+e_et[1,:])[mtotal*mrtot],facecolor = 'C4', alpha = 0.3)
plt.plot(ang_ind,et[mtotal*mint],c='C4',ls = '--')
plt.fill_between(ang_ind, (et-e_et[0,:])[mtotal*mint], (et+e_et[1,:])[mtotal*mint],facecolor = 'C4', alpha = 0.3)
plt.plot(ang_ind,et[mtotal*mext],c='C4',ls = '-.')
plt.fill_between(ang_ind, (et-e_et[0,:])[mtotal*mext], (et+e_et[1,:])[mtotal*mext],facecolor = 'C4', alpha = 0.3)
plt.xticks(ang_ind,angles)
plt.legend()
plt.savefig(folder+'plot_results/tot_inout_t.png')

plt.figure()
plt.plot(ang_ind,et[mmbin1*mrtot],c='C5',label = 'lM200_mean ='+str(lM200[mmbin1*mrtot][0]))
plt.fill_between(ang_ind, (et-e_et[0,:])[mmbin1*mrtot], (et+e_et[1,:])[mmbin1*mrtot],facecolor = 'C5', alpha = 0.3)
plt.plot(ang_ind,et[mmbin1*mint],c='C5',ls = '--')
plt.fill_between(ang_ind, (et-e_et[0,:])[mmbin1*mint], (et+e_et[1,:])[mmbin1*mint],facecolor = 'C5', alpha = 0.3)
plt.plot(ang_ind,et[mmbin1*mext],c='C5',ls = '-.')
plt.fill_between(ang_ind, (et-e_et[0,:])[mmbin1*mext], (et+e_et[1,:])[mmbin1*mext],facecolor = 'C5', alpha = 0.3)
plt.xticks(ang_ind,angles)
plt.legend()
plt.savefig(folder+'plot_results/LM_inout_t.png')

plt.figure()
plt.plot(ang_ind,et[mmbin2*mrtot],c='C6',label = 'lM200_mean ='+str(lM200[mmbin2*mrtot][0]))
plt.fill_between(ang_ind, (et-e_et[0,:])[mmbin2*mrtot], (et+e_et[1,:])[mmbin2*mrtot],facecolor = 'C6', alpha = 0.3)
plt.plot(ang_ind,et[mmbin2*mint],c='C6',ls = '--')
plt.fill_between(ang_ind, (et-e_et[0,:])[mmbin2*mint], (et+e_et[1,:])[mmbin2*mint],facecolor = 'C6', alpha = 0.3)
plt.plot(ang_ind,et[mmbin2*mext],c='C6',ls = '-.')
plt.fill_between(ang_ind, (et-e_et[0,:])[mmbin2*mext], (et+e_et[1,:])[mmbin2*mext],facecolor = 'C6', alpha = 0.3)
plt.xticks(ang_ind,angles)
plt.legend()
plt.savefig(folder+'plot_results/HM_inout_t.png')

plt.figure()
plt.plot(ang_ind,et[mtotal*mint],c='C4',label = 'lM200_mean ='+str(lM200[mtotal*mint][0]),ls = '--')
plt.fill_between(ang_ind, (et-e_et[0,:])[mtotal*mint], (et+e_et[1,:])[mtotal*mint],facecolor = 'C4', alpha = 0.3)
plt.plot(ang_ind,et[mmbin1*mint],c='C5',label = 'lM200_mean ='+str(lM200[mmbin1*mint][0]),ls = '--')
plt.fill_between(ang_ind, (et-e_et[0,:])[mmbin1*mint], (et+e_et[1,:])[mmbin1*mint],facecolor = 'C5', alpha = 0.3)
plt.plot(ang_ind,et[mmbin2*mint],c='C6',label = 'lM200_mean ='+str(lM200[mmbin2*mint][0]),ls = '--')
plt.fill_between(ang_ind, (et-e_et[0,:])[mmbin2*mint], (et+e_et[1,:])[mmbin2*mint],facecolor = 'C6', alpha = 0.3)
plt.xticks(ang_ind,angles)
plt.legend()
plt.savefig(folder+'plot_results/Mrin_t.png')

plt.figure()
plt.plot(ang_ind,et[mtotal*mext],c='C4',label = 'lM200_mean ='+str(lM200[mtotal*mext][0]),ls = '-.')
plt.fill_between(ang_ind, (et-e_et[0,:])[mtotal*mext], (et+e_et[1,:])[mtotal*mext],facecolor = 'C4', alpha = 0.3)
plt.plot(ang_ind,et[mmbin1*mext],c='C5',label = 'lM200_mean ='+str(lM200[mmbin1*mext][0]),ls = '-.')
plt.fill_between(ang_ind, (et-e_et[0,:])[mmbin1*mext], (et+e_et[1,:])[mmbin1*mint],facecolor = 'C5', alpha = 0.3)
plt.plot(ang_ind,et[mmbin2*mext],c='C6',label = 'lM200_mean ='+str(lM200[mmbin2*mext][0]),ls = '-.')
plt.fill_between(ang_ind, (et-e_et[0,:])[mmbin2*mext], (et+e_et[1,:])[mmbin2*mext],facecolor = 'C6', alpha = 0.3)
plt.xticks(ang_ind,angles)
plt.legend()
plt.savefig(folder+'plot_results/Mrout_t.png')

# REDSHIFT

plt.figure()
plt.plot(ang_ind,et[mtotal*mrtot],c='C4',label = 'z_mean ='+str(zmean[mtotal*mrtot][0]))
plt.fill_between(ang_ind, (et-e_et[0,:])[mtotal*mrtot], (et+e_et[1,:])[mtotal*mrtot],facecolor = 'C4', alpha = 0.3)
plt.plot(ang_ind,et[mmz1*mrtot],c='C5',label = 'z_mean ='+str(zmean[mmz1*mrtot][0]))
plt.fill_between(ang_ind, (et-e_et[0,:])[mmz1*mrtot], (et+e_et[1,:])[mmz1*mrtot],facecolor = 'C5', alpha = 0.3)
plt.plot(ang_ind,et[mmz2*mrtot],c='C6',label = 'z_mean ='+str(zmean[mmz2*mrtot][0]))
plt.fill_between(ang_ind, (et-e_et[0,:])[mmz2*mrtot], (et+e_et[1,:])[mmz2*mrtot],facecolor = 'C6', alpha = 0.3)
plt.xticks(ang_ind,angles)
plt.legend()
plt.savefig(folder+'plot_results/zrtot_t.png')

plt.figure()
plt.plot(ang_ind,et[mtotal*mrtot],c='C4',label = 'z_mean ='+str(zmean[mtotal*mrtot][0]))
plt.fill_between(ang_ind, (et-e_et[0,:])[mtotal*mrtot], (et+e_et[1,:])[mtotal*mrtot],facecolor = 'C4', alpha = 0.3)
plt.plot(ang_ind,et[mtotal*mint],c='C4',ls = '--')
plt.fill_between(ang_ind, (et-e_et[0,:])[mtotal*mint], (et+e_et[1,:])[mtotal*mint],facecolor = 'C4', alpha = 0.3)
plt.plot(ang_ind,et[mtotal*mext],c='C4',ls = '-.')
plt.fill_between(ang_ind, (et-e_et[0,:])[mtotal*mext], (et+e_et[1,:])[mtotal*mext],facecolor = 'C4', alpha = 0.3)
plt.xticks(ang_ind,angles)
plt.legend()

plt.figure()
plt.plot(ang_ind,et[mmz1*mrtot],c='C5',label = 'z_mean ='+str(zmean[mmz1*mrtot][0]))
plt.fill_between(ang_ind, (et-e_et[0,:])[mmz1*mrtot], (et+e_et[1,:])[mmz1*mrtot],facecolor = 'C5', alpha = 0.3)
plt.plot(ang_ind,et[mmz1*mint],c='C5',ls = '--')
plt.fill_between(ang_ind, (et-e_et[0,:])[mmz1*mint], (et+e_et[1,:])[mmz1*mint],facecolor = 'C5', alpha = 0.3)
plt.plot(ang_ind,et[mmz1*mext],c='C5',ls = '-.')
plt.fill_between(ang_ind, (et-e_et[0,:])[mmz1*mext], (et+e_et[1,:])[mmz1*mext],facecolor = 'C5', alpha = 0.3)
plt.xticks(ang_ind,angles)
plt.legend()
plt.savefig(folder+'plot_results/Lz_inout_t.png')

plt.figure()
plt.plot(ang_ind,et[mmz2*mrtot],c='C6',label = 'z_mean ='+str(zmean[mmz2*mrtot][0]))
plt.fill_between(ang_ind, (et-e_et[0,:])[mmz2*mrtot], (et+e_et[1,:])[mmz2*mrtot],facecolor = 'C6', alpha = 0.3)
plt.plot(ang_ind,et[mmz2*mint],c='C6',ls = '--')
plt.fill_between(ang_ind, (et-e_et[0,:])[mmz2*mint], (et+e_et[1,:])[mmz2*mint],facecolor = 'C6', alpha = 0.3)
plt.plot(ang_ind,et[mmz2*mext],c='C6',ls = '-.')
plt.fill_between(ang_ind, (et-e_et[0,:])[mmz2*mext], (et+e_et[1,:])[mmz2*mext],facecolor = 'C6', alpha = 0.3)
plt.xticks(ang_ind,angles)
plt.legend()
plt.savefig(folder+'plot_results/Hz_inout_t.png')

plt.figure()
plt.plot(ang_ind,et[mtotal*mint],c='C4',label = 'z_mean ='+str(zmean[mtotal*mint][0]),ls = '--')
plt.fill_between(ang_ind, (et-e_et[0,:])[mtotal*mint], (et+e_et[1,:])[mtotal*mint],facecolor = 'C4', alpha = 0.3)
plt.plot(ang_ind,et[mmz1*mint],c='C5',label = 'z_mean ='+str(zmean[mmz1*mint][0]),ls = '--')
plt.fill_between(ang_ind, (et-e_et[0,:])[mmz1*mint], (et+e_et[1,:])[mmz1*mint],facecolor = 'C5', alpha = 0.3)
plt.plot(ang_ind,et[mmz2*mint],c='C6',label = 'z_mean ='+str(zmean[mmz2*mint][0]),ls = '--')
plt.fill_between(ang_ind, (et-e_et[0,:])[mmz2*mint], (et+e_et[1,:])[mmz2*mint],facecolor = 'C6', alpha = 0.3)
plt.xticks(ang_ind,angles)
plt.legend()
plt.savefig(folder+'plot_results/zrin_t.png')

plt.figure()
plt.plot(ang_ind,et[mtotal*mext],c='C4',label = 'z_mean ='+str(zmean[mtotal*mext][0]),ls = '-.')
plt.fill_between(ang_ind, (et-e_et[0,:])[mtotal*mext], (et+e_et[1,:])[mtotal*mext],facecolor = 'C4', alpha = 0.3)
plt.plot(ang_ind,et[mmz1*mext],c='C5',label = 'z_mean ='+str(zmean[mmz1*mext][0]),ls = '-.')
plt.fill_between(ang_ind, (et-e_et[0,:])[mmz1*mext], (et+e_et[1,:])[mmz1*mext],facecolor = 'C5', alpha = 0.3)
plt.plot(ang_ind,et[mmz2*mext],c='C6',label = 'z_mean ='+str(zmean[mmz2*mext][0]),ls = '-.')
plt.fill_between(ang_ind, (et-e_et[0,:])[mmz2*mext], (et+e_et[1,:])[mmz2*mext],facecolor = 'C6', alpha = 0.3)
plt.xticks(ang_ind,angles)
plt.legend()
plt.savefig(folder+'plot_results/zrout_t.png')
