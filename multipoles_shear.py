import numpy as np
from pylab import *
from astropy.cosmology import LambdaCDM
from scipy.misc import derivative
from scipy import integrate
from profiles_fit import *
from multiprocessing import Pool
from multiprocessing import Process
import time
import sys, os

# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore
def enablePrint():
    sys.stdout = sys.__stdout__
#parameters

cvel = 299792458;   # Speed of light (m.s-1)
G    = 6.670e-11;   # Gravitational constant (m3.kg-1.s-2)
pc   = 3.085678e16; # 1 pc (m)
Msun = 1.989e30 # Solar mass (kg)

def multipole_shear(r,M200=1.e14,ellip=0.25,z=0.2,h=0.7,
					misscentred=False,s_off=0.4,
					components = ['t0','t','tcos','xcos'],
					verbose=True):

	'''
	Equations from van Uitert (vU, arXiv:1610.04226) for the 
	multipole expansion and from Ford et al. (F,2015) for
	the misscentring

	INPUTS:
	r               [Mpc] - Float or numpy_array
	                Distance to the centre of 
	                the gravitational potential
	
	OPTIONAL:
	M200            [M_sun] - Float/ Cluster Mass
	ellip           Float / Halo ellipticity defined as
	                (1-q)/(1+q) where q is the semi-axis
	                ratio (q < 1)
	z               Float/ Lens redshift
	zs              Float/ Source redshift
	h               Float/ Cosmological quantity
	misscentred     Float/ If True computes de misscentred quantities
					The misscentred is considered only in the 
					x - axis
	s_off           Float/ sigma_offset width of the distribution 
	                of cluster offsets (F_Eq11)
	components      List of misscentred components that are going
	                to be computed:
	                't0'   -> Gt0_off
	                't'    -> Gt_off
	                'tcos' -> Gt_off_cos
	                'xsin' -> Gx_off_sin
	verbose         Print time computing
	
	OUTPUT:
	output          Dictionary (the shear is scales with the
	                critical density Gamma_i = Sigma_cr*gamma_i
	                ---------- vU_Eq18 ---------
	                Gt0: Shear monopole
	                Gt2: Tangential component for the quadrupole
	                Gx2: Cross component for the quadrupole
	                --------- Misscentred terms --------
	                Gt0_off: (ellip = 0) F_Eq14
	                Gt_off: Integrated tangential component (vU_Eq17)
	                for considering kappa_offset
	                Gt_off_cos: Integrated tangential component times
	                cos(2theta) (vU_Eq17) for considering kappa_offset
	                Gx_off_sin: Integrated cross component times
	                sin(2theta) (vU_Eq17) for considering kappa_offset
	'''
	
	if not verbose:
		blockPrint()

	# Check if r is float or numpy array
	if not isinstance(r, (np.ndarray)):
		r = np.array([r])

	# Compute cosmological parameters
	cosmo = LambdaCDM(H0=h*100, Om0=0.3, Ode0=0.7)
	H        = cosmo.H(z).value/(1.0e3*pc) #H at z_pair s-1 
	roc      = (3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_pair (kg.m-3)
	roc_mpc  = roc*((pc*1.0e6)**3.0)
	
	# Compute R_200
	R200 = r200_nfw(M200,roc_mpc)

	# Scaling sigma_off
	s_off = s_off/h	
	
	#print '##################'
	#print '      CENTRED     '
	#print '##################'
	
	def Delta_Sigma(R):
		'''
		Density contraste for NFW
		
		'''
		
		#calculo de c usando la relacion de Duffy et al 2008
		
		M=((800.0*np.pi*roc_mpc*(R200**3))/(3.0*Msun))*h
		c=5.71*((M/2.e12)**-0.084)*((1.+z)**-0.47)
		# c = 5.0
		####################################################
		
		deltac=(200./3.)*( (c**3) / ( np.log(1.+c)- (c/(1+c)) ))
		x=np.round((R*c)/R200,12)
		m1= x< 1.0
		m2= x> 1.0 
		m3= (x == 1.0)
		
		try: 
			jota=np.zeros(len(x))
			atanh=np.arctanh(((1.0-x[m1])/(1.0+x[m1]))**0.5)
			jota[m1]=(4.0*atanh)/((x[m1]**2.0)*((1.0-x[m1]**2.0)**0.5)) \
				+ (2.0*np.log(x[m1]/2.0))/(x[m1]**2.0) - 1.0/(x[m1]**2.0-1.0) \
				+ (2.0*atanh)/((x[m1]**2.0-1.0)*((1.0-x[m1]**2.0)**0.5))    
			atan=np.arctan(((x[m2]-1.0)/(1.0+x[m2]))**0.5)
			jota[m2]=(4.0*atan)/((x[m2]**2.0)*((x[m2]**2.0-1.0)**0.5)) \
				+ (2.0*np.log(x[m2]/2.0))/(x[m2]**2.0) - 1.0/(x[m2]**2.0-1.0) \
				+ (2.0*atan)/((x[m2]**2.0-1.0)**1.5)
			jota[m3]=2.0*np.log(0.5)+5.0/3.0
		except:
			if m1:
				atanh=np.arctanh(((1.0-x[m1])/(1.0+x[m1]))**0.5)
				jota = (4.0*atanh)/((x[m1]**2.0)*((1.0-x[m1]**2.0)**0.5)) \
					+ (2.0*np.log(x[m1]/2.0))/(x[m1]**2.0) - 1.0/(x[m1]**2.0-1.0) \
					+ (2.0*atanh)/((x[m1]**2.0-1.0)*((1.0-x[m1]**2.0)**0.5))   
			if m2:		 
				atan=np.arctan(((x[m2]-1.0)/(1.0+x[m2]))**0.5)
				jota = (4.0*atan)/((x[m2]**2.0)*((x[m2]**2.0-1.0)**0.5)) \
					+ (2.0*np.log(x[m2]/2.0))/(x[m2]**2.0) - 1.0/(x[m2]**2.0-1.0) \
					+ (2.0*atan)/((x[m2]**2.0-1.0)**1.5)
			if m3:
				jota = 2.0*np.log(0.5)+5.0/3.0
	
	
			
		rs_m=(R200*1.e6*pc)/c
		kapak=((2.*rs_m*deltac*roc_mpc)*(pc**2/Msun))/((pc*1.0e6)**3.0)
		return kapak*jota
	
	
	def monopole(R):		
		'''
		Projected density for NFW
		
		'''		
		if not isinstance(R, (np.ndarray)):
			R = np.array([R])
		
		# m = R == 0.
		# R[m] = 1.e-8
		
		#calculo de c usando la relacion de Duffy et al 2008
		
		M=((800.0*np.pi*roc_mpc*(R200**3))/(3.0*Msun))*h
		c=5.71*((M/2.e12)**-0.084)*((1.+z)**-0.47)
		# c = 5.0
		####################################################
		
		deltac=(200./3.)*( (c**3) / ( np.log(1.+c)- (c/(1+c)) ))

		x=(R*c)/R200
		m1 = x <= (1.0-1.e-4)
		m2 = x >= (1.0+1.e-4)
		m3 = (x == 1.0)
		m4 = (~m1)*(~m2)*(~m3)
	
		jota  = np.zeros(len(x))
		atanh = np.arctanh(np.sqrt((1.0-x[m1])/(1.0+x[m1])))
		jota[m1] = (1./(x[m1]**2-1.))*(1.-(2./np.sqrt(1.-x[m1]**2))*atanh) 
	
		atan = np.arctan(((x[m2]-1.0)/(1.0+x[m2]))**0.5)
		jota[m2] = (1./(x[m2]**2-1.))*(1.-(2./np.sqrt(x[m2]**2 - 1.))*atan) 
	
		jota[m3] = 1./3.
		
		x1 = 1.-1.e-4
		atanh1 = np.arctanh(np.sqrt((1.0-x1)/(1.0+x1)))
		j1 = (1./(x1**2-1.))*(1.-(2./np.sqrt(1.-x1**2))*atanh1) 
		
		x2 = 1.+1.e-4
		atan2 = np.arctan(((x2-1.0)/(1.0+x2))**0.5)
		j2 = (1./(x2**2-1.))*(1.-(2./np.sqrt(x2**2 - 1.))*atan2) 
		
		jota[m4] = np.interp(x[m4].astype(float64),[x1,x2],[j1,j2])
					
		rs_m=(R200*1.e6*pc)/c
		kapak=((2.*rs_m*deltac*roc_mpc)*(pc**2/Msun))/((pc*1.0e6)**3.0)
		return kapak*jota

	def quadrupole(R):
		'''
		Quadrupole term defined as (d(Sigma)/dr)*r
		
		'''		
		m0p = derivative(monopole,R,dx=1e-5)
		return m0p*R

	def psi2(R):
		'''
		vU_Eq10
		
		'''
		argumento = lambda x: (x**3)*monopole(x)
		integral = integrate.quad(argumento, 0, R)[0]
		return integral*(-2./(R**2))
	
	vecpsi2 = np.vectorize(psi2)
		
	#print '##################'
	#print '    MISCENTRED    '
	#print '##################'

	def P_Roff(Roff):
		'''
		F_Eq11
		
		'''		
		return abs((Roff/s_off**2)*np.exp(-0.5*(Roff/s_off)**2))*0.5
	
	def monopole_off(R,theta):
		'''
		F_Eq12
		
		'''				
		def moff(x):
			return monopole(R**2+x**2-2*x*R*np.cos(theta))*P_Roff(x)
		argumento = lambda x: moff(x)
		integral1  = integrate.quad(argumento, -1.*np.inf, 0)[0]
		integral2  = integrate.quad(argumento, 0., R)[0]
		integral3  = integrate.quad(argumento, R, np.inf)[0]
		return integral1 + integral2 + integral3
	vec_moff = np.vectorize(monopole_off)

	def Delta_Sigma_off(R,theta):
		'''
		F_Eq14
		
		'''						
		argumento = lambda x: monopole_off(x,theta)*x
		integral  = integrate.quad(argumento, 0, R)[0]
		DS_off    = (2./R**2)*integral - monopole_off(R,theta)
		return DS_off

	vec_DSoff = np.vectorize(Delta_Sigma_off)
	
	def quadrupole_off(R,theta):

		def rp(roff):
			return R**2+roff**2-2*roff*R*np.cos(theta)
			
		def q_off(roff):
			
			if rp(roff) > 1.e-4:
				return quadrupole(rp(roff))*P_Roff(roff)
			else:
				# for larger precission interpolates
				b = -2*R*np.cos(theta)
				c = R**2-1.e-4
				R1 = np.round(((-b - np.sqrt(b**2 - 4.*c))*0.5),6)
				R2 = np.round(((-b + np.sqrt(b**2 - 4.*c))*0.5),6)
				R1 = min(R1,R2)		
				R2 = max(R1,R2)

				q1 = (quadrupole(rp(R1))*P_Roff(R1))[0]
				q2 = (quadrupole(rp(R2))*P_Roff(R2))[0]
				
				return np.interp(roff,[R1,R2],[q1,q2])

		argumento = lambda x: q_off(x)
		integral10  = integrate.quad(argumento, -1.*np.inf, 0)[0]
		integral20  = integrate.quad(argumento, 0., R)[0]
		integral30  = integrate.quad(argumento, R, np.inf)[0]

		return integral10 + integral20 + integral30	
	vec_qoff = np.vectorize(quadrupole_off)
	
	def psi2_off(R,theta):
		def arg(x):
			return (x**3)*monopole_off(x,theta)
		argumento = lambda x: arg(x)
		integral = integrate.quad(argumento, 0, R)[0]
		return integral*(-2./(R**2))
		
	vecpsi2_off = np.vectorize(psi2_off)	

	def quantities_centred(r):
		
		gt0 = Delta_Sigma(r)
		monopole_r = monopole(r)
		quadrupole_r = quadrupole(r)
		print 'computing psi2 centred'	
		psi2_r = vecpsi2(r)
				
		return gt0,monopole_r,quadrupole_r,psi2_r
		
	def quantities_misscentred(r):
		print 'computing misscentred profile'
		gamma_t0_off = []
		gamma_t_off0 = []
		gamma_t_off = []
		gamma_x_off = []
		for R in r:
			
			if 't0' in components:
				print 'computing DS_t0_off'
				t1 = time.time()
				def DS_t_off0(theta):
					gamma_t0 = Delta_Sigma_off(R,theta)
					return gamma_t0
				argumento = lambda x: DS_t_off0(x)
				integral  = integrate.quad(argumento, 0, 2.*np.pi,points=[np.pi], epsabs=1.e-04, epsrel=1.e-04)[0]
				gamma_t0_off = np.append(gamma_t0_off,integral/(2.*np.pi))
				t2 = time.time()
				print (t2-t1)/60.			

			def DS_t_off(theta):
				gamma_t0 = Delta_Sigma_off(R,theta)
				gamma_t2 = ((-6*psi2_off(R,theta)/R**2) 
				            - 2.*monopole_off(R,theta) 
				            + quadrupole_off(R,theta))
				return gamma_t0 + ellip*gamma_t2*np.cos(2.*theta)

			if 't' in components:
				print 'computing DS_t_off'
				print R,M200
				t1 = time.time()
				argumento = lambda x: DS_t_off(x)
				integral  = integrate.quad(argumento, 0., 2.*np.pi,points=[np.pi], epsabs=1.e-04, epsrel=1.e-04)[0]
				print integral
				gamma_t_off0 = np.append(gamma_t_off0,integral/(2.*np.pi))
				t2 = time.time()
			 	print (t2-t1)/60.
			
			if 'tcos' in components:	
				print 'computing DS_t_off_cos'
				t1 = time.time()
				argumento = lambda x: DS_t_off(x)*np.cos(2.*x)
				integral  = integrate.quad(argumento, 0., 2.*np.pi,points=[np.pi], epsabs=1.e-04, epsrel=1.e-04)[0]
				gamma_t_off = np.append(gamma_t_off,integral/np.pi)
				t2 = time.time()
			 	print (t2-t1)/60.
			
			if 'xcos' in components:
				print 'computing DS_x_off_sin'
				t1 = time.time()
				def DS_x_off(theta):
					gamma_x2 = ((-6*psi2_off(R,theta)/R**2) 
								- 4.*monopole_off(R,theta))
					return ellip*gamma_x2*np.sin(2.*theta)
				argumento = lambda x: DS_x_off(x)*np.sin(2.*x)
				integral  = integrate.quad(argumento, 0, 2.*np.pi,points=[np.pi], epsabs=1.e-04, epsrel=1.e-04)[0]
				gamma_x_off = np.append(gamma_x_off,integral/np.pi)	
				t2 = time.time()
				print (t2-t1)/60.
					
		return gamma_t0_off, gamma_t_off0, gamma_t_off, gamma_x_off
		
	
	gt0,m,q,p2 = quantities_centred(r)
	
	gt2 = ellip*((-6*p2/r**2) - 2.*m + q)
	gx2 = ellip*((-6*p2/r**2) - 4.*m)
	
	output = {'Gt0':gt0,'Gt2':gt2,'Gx2':gx2}
	
	if misscentred:
		gt0_off, gt_off0, gt_off, gx_off = quantities_misscentred(r)	
		output.update({'Gt0_off':gt0_off,'Gt_off':gt_off0,'Gt_off_cos':gt_off,'Gx_off_sin':gx_off})

	enablePrint()
		
	return output

# '''

def multipole_shear_unpack(minput):
	return multipole_shear(*minput)
	
def multipole_shear_parallel(r,M200=1.e14,ellip=0.25,z=0.2,
							 h=0.7,misscentred=False,
							 s_off=0.4,components = ['t0','t','tcos','xcos'],
							 verbose = True, ncores=2):
	
	if ncores > len(r):
		ncores = len(r)
	
	
	slicer = int(round(len(r)/float(ncores), 0))
	slices = ((np.arange(ncores-1)+1)*slicer).astype(int)
	slices = slices[(slices <= len(r))]
	r_splitted = np.split(r,slices)
	
	ncores = len(r_splitted)
	
	M200  = np.ones(ncores)*M200
	ellip = np.ones(ncores)*ellip
	z     = np.ones(ncores)*z
	h     = np.ones(ncores)*h
	miss  = np.ones(ncores,dtype=bool)*misscentred
	s_off = np.ones(ncores)*s_off
	comp  = [components]*ncores
	v  = np.ones(ncores,dtype=bool)*verbose
	
	entrada = np.array([r_splitted,M200,ellip,z,h,miss,s_off,comp,v]).T
	
	pool = Pool(processes=(ncores))
	salida=np.array(pool.map(multipole_shear_unpack, entrada))
	pool.terminate()

	gt0 = []
	gt2 = []
	gx2 = []
	
	if misscentred:
		gt0_off = []
		gt_off0 = []
		gt_off  = []
		gx_off  = []

	for s in salida:
		gt0 = np.append(gt0,s['Gt0'])
		gt2 = np.append(gt2,s['Gt2'])
		gx2 = np.append(gx2,s['Gx2'])
		if misscentred:
			gt0_off = np.append(gt0_off,s['Gt0_off'])
			gt_off0 = np.append(gt_off0,s['Gt_off'])
			gt_off  = np.append(gt_off,s['Gt_off_cos'])
			gx_off  = np.append(gx_off,s['Gx_off_sin'])
			
	output = {'Gt0':gt0,'Gt2':gt2,'Gx2':gx2}
	
	if misscentred:
		output.update({'Gt0_off':gt0_off,'Gt_off':gt_off0,'Gt_off_cos':gt_off,'Gx_off_sin':gx_off})

	return output


def model_Gamma(multipole_out,component = 't', misscentred = False, pcc = 0.8):
	
	if misscentred:
		if component == 't':
			G = pcc*multipole_out['Gt0'] + (1-pcc)*multipole_out['Gt_off']
		if component == 'tcos':
			G = pcc*multipole_out['Gt2'] + (1-pcc)*multipole_out['Gt_off_cos']
		if component == 'xsin':
			G = pcc*multipole_out['Gx2'] + (1-pcc)*multipole_out['Gx_off_sin']
	else:
		if component == 't':
			G = multipole_out['Gt0'] 
		if component == 'tcos':
			G = multipole_out['Gt2'] 
		if component == 'xsin':
			G = multipole_out['Gx2'] 
			
	return G

'''
def multipole_clampitt(r,M200=1.e14,z=0.2,zs=0.35,
					   h=0.7,misscentred=False,s_off=0.4):

	cosmo = LambdaCDM(H0=h*100, Om0=0.3, Ode0=0.7)
	H        = cosmo.H(z).value/(1.0e3*pc) #H at z_pair s-1 
	roc      = (3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_pair (kg.m-3)
	roc_mpc  = roc*((pc*1.0e6)**3.0)
	
	
	R200 = r200_nfw(M200,roc_mpc)
	
	s_off = s_off/h	
	
	############  COMPUTING S_crit
	
	Dl    = cosmo.angular_diameter_distance(z).value*1.e6*pc
	Ds    = cosmo.angular_diameter_distance(zs).value*1.e6*pc
	Dls   = cosmo.angular_diameter_distance_z1z2(z,zs).value*1.e6*pc
	
	Sc = ((((cvel**2.0)/(4.0*np.pi*G*Dl))*(1./(Dls/Ds)))*(pc**2/Msun))
	sigma_c = np.zeros(len(r))
	sigma_c.fill(Sc)


	#print '##################'
	#print '      CENTRED     '
	#print '##################'
	
	def Delta_Sigma(R):
		
		#calculo de c usando la relacion de Duffy et al 2008
		
		M=((800.0*np.pi*roc_mpc*(R200**3))/(3.0*Msun))*h
		c=5.71*((M/2.e12)**-0.084)*((1.+z)**-0.47)
		# c = 5.0
		####################################################
		
		deltac=(200./3.)*( (c**3) / ( np.log(1.+c)- (c/(1+c)) ))
		x=(R*c)/R200
		m1= x< 1.0
		m2= x> 1.0 
		m3= (x == 1.0)
		
		try: 
			jota=np.zeros(len(x))
			atanh=np.arctanh(((1.0-x[m1])/(1.0+x[m1]))**0.5)
			jota[m1]=(4.0*atanh)/((x[m1]**2.0)*((1.0-x[m1]**2.0)**0.5)) \
				+ (2.0*np.log(x[m1]/2.0))/(x[m1]**2.0) - 1.0/(x[m1]**2.0-1.0) \
				+ (2.0*atanh)/((x[m1]**2.0-1.0)*((1.0-x[m1]**2.0)**0.5))    
			atan=np.arctan(((x[m2]-1.0)/(1.0+x[m2]))**0.5)
			jota[m2]=(4.0*atan)/((x[m2]**2.0)*((x[m2]**2.0-1.0)**0.5)) \
				+ (2.0*np.log(x[m2]/2.0))/(x[m2]**2.0) - 1.0/(x[m2]**2.0-1.0) \
				+ (2.0*atan)/((x[m2]**2.0-1.0)**1.5)
			jota[m3]=2.0*np.log(0.5)+5.0/3.0
		except:
			if m1:
				atanh=np.arctanh(((1.0-x[m1])/(1.0+x[m1]))**0.5)
				jota = (4.0*atanh)/((x[m1]**2.0)*((1.0-x[m1]**2.0)**0.5)) \
					+ (2.0*np.log(x[m1]/2.0))/(x[m1]**2.0) - 1.0/(x[m1]**2.0-1.0) \
					+ (2.0*atanh)/((x[m1]**2.0-1.0)*((1.0-x[m1]**2.0)**0.5))   
			if m2:		 
				atan=np.arctan(((x[m2]-1.0)/(1.0+x[m2]))**0.5)
				jota = (4.0*atan)/((x[m2]**2.0)*((x[m2]**2.0-1.0)**0.5)) \
					+ (2.0*np.log(x[m2]/2.0))/(x[m2]**2.0) - 1.0/(x[m2]**2.0-1.0) \
					+ (2.0*atan)/((x[m2]**2.0-1.0)**1.5)
			if m3:
				jota = 2.0*np.log(0.5)+5.0/3.0
	
	
			
		rs_m=(R200*1.e6*pc)/c
		kapak=((2.*rs_m*deltac*roc_mpc)*(pc**2/Msun))/((pc*1.0e6)**3.0)
		return kapak*jota
	
	
	def monopole(R):
		
		if not isinstance(R, (np.ndarray)):
			R = np.array([R])
		
		# m = R == 0.
		# R[m] = 1.e-8
		
		#calculo de c usando la relacion de Duffy et al 2008
		
		M=((800.0*np.pi*roc_mpc*(R200**3))/(3.0*Msun))*h
		c=5.71*((M/2.e12)**-0.084)*((1.+z)**-0.47)
		# c = 5.0
		####################################################
		
		deltac=(200./3.)*( (c**3) / ( np.log(1.+c)- (c/(1+c)) ))
		x=(R*c)/R200
		m1= x< 1.0
		m2= x> 1.0 
		m3= (x == 1.0)
	
		jota  = np.zeros(len(x))
		atanh = np.arctanh(np.sqrt((1.0-x[m1])/(1.0+x[m1])))
		jota[m1] = (1./(x[m1]**2-1.))*(1.-(2./np.sqrt(1.-x[m1]**2))*atanh) 
	
		atan = np.arctan(((x[m2]-1.0)/(1.0+x[m2]))**0.5)
		jota[m2] = (1./(x[m2]**2-1.))*(1.-(2./np.sqrt(x[m2]**2 - 1.))*atan) 
	
		jota[m3] = 1./3.
					
		rs_m=(R200*1.e6*pc)/c
		kapak=((2.*rs_m*deltac*roc_mpc)*(pc**2/Msun))/((pc*1.0e6)**3.0)
		return kapak*jota

	def quadrupole(R):
		m0p = derivative(monopole,R,dx=1e-6)
		return m0p*R
	
	def I1(R):
		argumento = lambda x: (x**3)*quadrupole(x)
		integral = integrate.quad(argumento, 0, R)[0]
		return integral*(3./(R**4))
	
	def I2(R):	
		argumento = lambda x: (quadrupole(x)/x)
		integral = integrate.quad(argumento, R, np.inf)[0]
		return integral
	
	vecI1 = np.vectorize(I1)
	vecI2 = np.vectorize(I2)
	
		
	#print '##################'
	#print '    MISCENTRED    '
	#print '##################'

	def P_Roff(Roff):
		return abs((Roff/s_off**2)*np.exp(-0.5*(Roff/s_off)**2))
	
	def monopole_off(R,theta):
		argumento = lambda x: monopole(R**2+x**2-2*x*R*np.cos(theta))*P_Roff(x)
		integral1  = integrate.quad(argumento, -1.*np.inf, R)[0]
		integral2  = integrate.quad(argumento, 0., R)[0]
		integral3  = integrate.quad(argumento, R, np.inf)[0]
		return integral1 + integral2 + integral3
	vec_moff = np.vectorize(monopole_off)
	
	def Delta_Sigma_off(R,theta):
		argumento = lambda x: monopole_off(x,theta)*x
		integral  = integrate.quad(argumento, 0, R)[0]
		DS_off    = (2./R**2)*integral - monopole_off(x,theta)
	vec_DSoff = np.vectorize(Delta_Sigma_off)
	
	def quadrupole_off(R,theta):
		def roff(x,R,theta):
			return np.round(R**2+x**2-2*x*R*np.cos(theta),6)
		def qoff(x,R,theta):
			return quadrupole(roff(x,R,theta))*P_Roff(x)
		argumento = lambda x: qoff(x,R,theta)
		integral1  = integrate.quad(argumento, -1.*np.inf, 0)[0]
		integral2  = integrate.quad(argumento, 0., R)[0]
		integral3  = integrate.quad(argumento, R, np.inf)[0]
		return integral1 + integral2 + integral3	
	vec_qoff = np.vectorize(quadrupole_off)
	
	
	def I1_off(R,theta):
		if theta != 0. and theta != np.pi:
			argumento = lambda x: (x**3)*quadrupole_off(x,theta)
			integral = integrate.quad(argumento, 0, R)[0]
			return integral*(3./(R**4))
		else:
			return 0.
		
	def I2_off(R,theta):	
		if theta != 0. and theta != np.pi:
			argumento = lambda x: quadrupole_off(x,theta)/x
			integral = integrate.quad(argumento, R, np.inf)[0]
			return integral
		else:
			return 0.
	vecI1_off = np.vectorize(I1_off)		
	vecI2_off = np.vectorize(I2_off)
			

	def quantities_centred(r):
		
		# optimize using unique r
		r,c = np.unique(r,return_counts=True)
		
		monopole_r = monopole(r)
		quadrupole_r = quadrupole(r)
		print 'computing I1 centred'	
		I1r = vecI1(r)
		print 'computing I2 centred'	
		I2r = vecI2(r)
		
		monopole_r   = np.repeat(monopole_r,c)
		quadrupole_r = np.repeat(quadrupole_r,c)
		I1r          = np.repeat(I1r,c)
		I2r          = np.repeat(I2r,c)
		
		return monopole_r,quadrupole_r,I1r,I2r
		
	def quantities_misscentred(r,theta):
		print 'computing monopole misscentred'	
		monopole_off_r = vec_moff(r,theta) 
		print 'computing quadrupole misscentred'	
		quadrupole_off_r = vec_qoff(r,theta)
		print 'computing I1 misscentred'	
		I1r_off = vecI1_off(r,theta)
		print 'computing I2 misscentred'	
		I2r_off = vecI2_off(r,theta)
		return monopole_off_r,quadrupole_off_r,I1r_off,I2r_off
	
	if misscentred:
		m,q,I1r,I2r = quantities_misscentred(r,theta)
	else:
		m,q,I1r,I2r = quantities_centred(r)
	
	output = {'monopole':m,'quadrupole':q,'I1':I1r,'I2':I2r}
		
	return output
'''
