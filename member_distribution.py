import numpy as np
from profiles_fit import SIGMA_nfw
from profiles_fit import r200_nfw
from astropy.cosmology import LambdaCDM

def momentos(dx,dy,w):
     
     Q11  = np.sum((dx**2)*w)/np.sum(w)
     Q22  = np.sum((dy**2)*w)/np.sum(w)
     Q12  = np.sum((dx*dy)*w)/np.sum(w)
     E1 = (Q11-Q22)/(Q11+Q22)
     E2 = (2.*Q12)/(Q11+Q22)
     e = np.sqrt(E1**2 + E2**2)
     theta = np.arctan2(E2,E1)/2.
     return e,theta
     


def ellip_fit(e,theta):
     
     
     a = np.linspace(-10,10,50)
     x,y = np.meshgrid(a,a)
     
     a = 5.
     b = a*((1. - e)/(1. + e))
     
     print (a-b)/(a+b)
     t = np.deg2rad(-1.*theta)
     
     mask = ((x**2/a**2) + (y**2/b**2)) < 1.
     
     x_rot = (x*np.cos(t) + y*np.sin(t)) 
     y_rot = (-1.*x*np.sin(t) + y*np.cos(t))
     
     w = np.ones(len(x_rot[mask]))
     e,ang=momentos(x_rot[mask],y_rot[mask],w)
     
     
     return e,np.rad2deg(ang)
     
def NFW_fit(e,theta,Lambda,z,h=0.7):
     
     
     a = np.linspace(-10,10,100)
     x,y = np.meshgrid(a,a)
     
     q = (1. - e)/(1. + e)
     
     r = np.sqrt(q*(x**2) + (y**2)/q)     
     
     # SIMET RELATION FOR M200
     M0    = (2.21e14)/h
     alpha = 1.33
     M200 = M0*((Lambda/40.)**alpha)
     
     # Compute cosmological parameters
     cosmo = LambdaCDM(H0=h*100, Om0=0.3, Ode0=0.7)
     H        = cosmo.H(z).value/(1.0e3*pc) #H at z_pair s-1 
     roc      = (3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_pair (kg.m-3)
     roc_mpc  = roc*((pc*1.0e6)**3.0)
     
     # Compute R_200
     R200 = r200_nfw(M200,roc_mpc)
     
     # S_nfw R_200
     S_nfw = SIGMA_nfw((r,sigma_c,z,h),R200)
     
     
     x_rot = (x*np.cos(t) + y*np.sin(t)) 
     y_rot = (-1.*x*np.sin(t) + y*np.cos(t))
     
     w = np.ones(len(x_rot[mask]))
     e,ang=momentos(x_rot[mask],y_rot[mask],w)
     
     
     return e,np.rad2deg(ang)
