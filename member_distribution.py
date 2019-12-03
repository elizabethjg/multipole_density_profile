import numpy as np
# from pylab import *
from astropy.io import fits
#parameters
from astropy.cosmology import LambdaCDM
from astropy.wcs import WCS
cosmo = LambdaCDM(H0=70., Om0=0.3, Ode0=0.7)

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
     

