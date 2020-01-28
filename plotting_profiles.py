import numpy as np
from pylab import *

f, ax1 = plt.subplots()
f, ax2 = plt.subplots()
f, ax3 = plt.subplots()


for q in np.arange(0.1,1.,0.2):

    out = np.loadtxt('multipoles_'+str('%2.1f' % q)+'.out')
    
    r = out[:,0]
    gt0 = out[:,1]
    gt2 = out[:,2]
    gx2 = out[:,3]
    
    gt0_off = out[:,4]
    gt_off  = out[:,5]
    gt_cos_off = out[:,6]
    gx_sin_off = out[:,7]
    
    pcc = 0.8
    
    gt = pcc*gt0 + (1-pcc)*gt0_off
    gt_e = pcc*gt0 + (1-pcc)*gt_off
    

    ax1.set_title('Gamma_t')
    ax1.plot(r,gt0,'C0',label = 'Only monopole',alpha = 1.-(q/2.))
    ax1.plot(r,gt,'C2',label = 'Monopole + misscentred (e=0)',alpha = 1.-(q/2.))
    ax1.plot(r,gt_e,'C1',label = 'Monopole + misscentred (e=0.25)',alpha = 1.-(q/2.))
    ax1.plot(r,(1-pcc)*gt_off,'C1--',label = 'Monopole + misscentred (e=0.25)',alpha = 1.-(q/2.))
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel(r'$\Gamma_t0$')
    ax1.set_xlabel('r [Mpc]')
    ax1.axis([0,1.5,1,350.])    
    # ax1.legend()

    
    # FOR ellip = 0.25
    
    gt_cos = pcc*gt2 + (1-pcc)*gt_cos_off
    

    ax2.set_title('Gamma_t_cos(2theta)',alpha = 1.-(q/2.))
    ax2.plot(r,gt2,'C0',label = 'Only quadrupole',alpha = 1.-(q/2.))
    ax2.plot(r,gt_cos,'C1',label = 'quadrupole + misscentred',alpha = 1.-(q/2.))
    ax2.plot(r,(1-pcc)*gt_cos_off,'C1--',label = 'quadrupole + misscentred',alpha = 1.-(q/2.))
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel(r'$\Gamma_t2$')
    ax2.set_xlabel('r [Mpc]')
    ax2.axis([0,1.5,1,350.])
    # ax2.legend()
    
    
    gx_sin = pcc*gx2 + (1-pcc)*gx_sin_off
    

    ax3.set_title('Gamma_x_sin(2theta)',alpha = 1.-(q/2.))
    ax3.plot(r,gx2,'C0',label = 'Only quadrupole',alpha = 1.-(q/2.))
    ax3.plot(r,gx_sin,'C1',label = 'quadrupole + misscentred',alpha = 1.-(q/2.))
    ax3.plot(r,(1-pcc)*gx_sin_off,'C1--',label = 'quadrupole + misscentred',alpha = 1.-(q/2.))
    ax3.set_xscale('log')
    ax3.set_ylabel(r'$\Gamma_x2$')
    ax3.set_xlabel('r [Mpc]')
    ax3.axis([0,1.5,-50.,30.])
    # ax3.legend()

f, ax1 = plt.subplots()
f, ax2 = plt.subplots()
f, ax3 = plt.subplots()
q = 0.

for logM in np.arange(12.5,15.5,0.5):

    out = np.loadtxt('multipoles_'+str(logM)+'.out')
    
    r = out[:,0]
    gt0 = out[:,1]
    gt2 = out[:,2]
    gx2 = out[:,3]
    
    gt0_off = out[:,4]
    gt_off  = out[:,5]
    gt_cos_off = out[:,6]
    gx_sin_off = out[:,7]
    
    pcc = 0.8
    
    gt = pcc*gt0 + (1-pcc)*gt0_off
    gt_e = pcc*gt0 + (1-pcc)*gt_off
    

    ax1.set_title('Gamma_t')
    ax1.plot(r,gt0,'C0',label = 'Only monopole',alpha = 1.-(q/2.))
    ax1.plot(r,gt,'C2',label = 'Monopole + misscentred (e=0)',alpha = 1.-(q/2.))
    ax1.plot(r,gt_e,'C1',label = 'Monopole + misscentred (e=0.25)',alpha = 1.-(q/2.))
    ax1.plot(r,(1-pcc)*gt_off,'C1--',label = 'Monopole + misscentred (e=0.25)',alpha = 1.-(q/2.))
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel(r'$\Gamma_t0$')
    ax1.set_xlabel('r [Mpc]')
    ax1.axis([0,1.5,1,350.])    
    # ax1.legend()

    
    # FOR ellip = 0.25
    
    gt_cos = pcc*gt2 + (1-pcc)*gt_cos_off
    

    ax2.set_title('Gamma_t_cos(2theta)',alpha = 1.-(q/2.))
    ax2.plot(r,gt2,'C0',label = 'Only quadrupole',alpha = 1.-(q/2.))
    ax2.plot(r,gt_cos,'C1',label = 'quadrupole + misscentred',alpha = 1.-(q/2.))
    ax2.plot(r,(1-pcc)*gt_cos_off,'C1--',label = 'quadrupole + misscentred',alpha = 1.-(q/2.))
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel(r'$\Gamma_t2$')
    ax2.set_xlabel('r [Mpc]')
    ax2.axis([0,1.5,1,350.])
    # ax2.legend()
    
    
    gx_sin = pcc*gx2 + (1-pcc)*gx_sin_off
    

    ax3.set_title('Gamma_x_sin(2theta)',alpha = 1.-(q/2.))
    ax3.plot(r,gx2,'C0',label = 'Only quadrupole',alpha = 1.-(q/2.))
    ax3.plot(r,gx_sin,'C1',label = 'quadrupole + misscentred',alpha = 1.-(q/2.))
    ax3.plot(r,(1-pcc)*gx_sin_off,'C1--',label = 'quadrupole + misscentred',alpha = 1.-(q/2.))
    ax3.set_xscale('log')
    ax3.set_ylabel(r'$\Gamma_x2$')
    ax3.set_xlabel('r [Mpc]')
    ax3.axis([0,1.5,-50.,30.])

print 'END OF PROGRAM'
