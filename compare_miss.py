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

out      = np.loadtxt(folder+'table_mcmc_centred.out',dtype='str').T
out_miss = np.loadtxt(folder+'table_mcmc_misscentred.out',dtype='str').T

samples = out[0]
M200    = out[1].astype(float)
lM200   = np.log10(M200*1.e14)
zmean   = out[2].astype(float)
rin     = out[3].astype(float)
rout    = out[4].astype(float)
eb      = out[5].astype(float)
e_eb    = np.array([out[6].astype(float),out[7].astype(float)])

eb_miss      = out_miss[5].astype(float)
e_eb_miss    = np.array([out_miss[6].astype(float),out_miss[7].astype(float)])


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
plt.plot(eb[mtotal*mrtot],eb_miss[mtotal*mrtot],'C0.',label='total')
plt.plot(eb[mtotal*mint],eb_miss[mtotal*mint],'C0x')
plt.plot(eb[mtotal*mext],eb_miss[mtotal*mext],'C0D')
plt.errorbar(eb[mtotal],eb_miss[mtotal],xerr=e_eb[:,mtotal],yerr=e_eb_miss[:,mtotal],fmt='none',ecolor='C0')
plt.plot([0,0.4],[0,0.4],'k--')
plt.legend()
plt.savefig(folder+'plot_results_misscentred/comparison_total.png')

plt.figure()
plt.plot(eb[mmbin1*mrtot],eb_miss[mmbin1*mrtot],'C1.',label='bin1')
plt.plot(eb[mmbin1*mint],eb_miss[mmbin1*mint],'C1x')
plt.plot(eb[mmbin1*mext],eb_miss[mmbin1*mext],'C1D')
plt.errorbar(eb[mmbin1],eb_miss[mmbin1],xerr=e_eb[:,mmbin1],yerr=e_eb_miss[:,mmbin1],fmt='none',ecolor='C1')
plt.plot([0,0.4],[0,0.4],'k--')
plt.legend()
plt.savefig(folder+'plot_results_misscentred/comparison_bin1.png')

plt.figure()
plt.plot(eb[mmbin2*mrtot],eb_miss[mmbin2*mrtot],'C2.',label='bin2')
plt.plot(eb[mmbin2*mint],eb_miss[mmbin2*mint],'C2x')
plt.plot(eb[mmbin2*mext],eb_miss[mmbin2*mext],'C2D')
plt.errorbar(eb[mmbin2],eb_miss[mmbin2],xerr=e_eb[:,mmbin2],yerr=e_eb_miss[:,mmbin2],fmt='none',ecolor='C2')
plt.plot([0,0.4],[0,0.4],'k--')
plt.legend()
plt.savefig(folder+'plot_results_misscentred/comparison_bin2.png')

plt.figure()
plt.plot(eb[mmz1*mrtot],eb_miss[mmz1*mrtot],'C3.',label='z1')
plt.plot(eb[mmz1*mint],eb_miss[mmz1*mint],'C3x')
plt.plot(eb[mmz1*mext],eb_miss[mmz1*mext],'C3D')
plt.errorbar(eb[mmz1],eb_miss[mmz1],xerr=e_eb[:,mmz1],yerr=e_eb_miss[:,mmz1],fmt='none',ecolor='C3')
plt.plot([0,0.4],[0,0.4],'k--')
plt.legend()
plt.savefig(folder+'plot_results_misscentred/comparison_z1.png')

plt.figure()
plt.plot(eb[mmz2*mrtot],eb_miss[mmz2*mrtot],'C4.',label='z2')
plt.plot(eb[mmz2*mint],eb_miss[mmz2*mint],'C4x')
plt.plot(eb[mmz2*mext],eb_miss[mmz2*mext],'C4D')
plt.errorbar(eb[mmz2],eb_miss[mmz2],xerr=e_eb[:,mmz2],yerr=e_eb_miss[:,mmz2],fmt='none',ecolor='C4')
plt.plot([0,0.4],[0,0.4],'k--')
plt.legend()
plt.savefig(folder+'plot_results_misscentred/comparison_z2.png')
