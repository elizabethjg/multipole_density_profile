import numpy as np
from pylab import *
from astropy.cosmology import LambdaCDM
from scipy.misc import derivative
from scipy import integrate
from profiles_fit import *
from multiprocessing import Pool
from multiprocessing import Process

#parameters

cvel = 299792458;   # Speed of light (m.s-1)
G    = 6.670e-11;   # Gravitational constant (m3.kg-1.s-2)
pc   = 3.085678e16; # 1 pc (m)
Msun = 1.989e30 # Solar mass (kg)

def multipole_clampitt(r,theta,M200=1.e14,z=0.2,zs=0.35,
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
		argumento = lambda x: quadrupole(R**2+x**2-2*x*R*np.cos(theta))*P_Roff(x)
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


def multipole_vanUitert(r,theta,M200=1.e14,z=0.2,zs=0.35,
					   h=0.7,misscentred=False,s_off=0.4,ncores=20):

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
		m0p = derivative(monopole,R,dx=1e-12)
		return m0p*R

	def psi2(R):
		argumento = lambda x: (x**3)*monopole(x)
		integral = integrate.quad(argumento, 0, R)[0]
		return integral*(-2./(R**2))
		
	vecpsi2 = np.vectorize(psi2)
		
	#print '##################'
	#print '    MISCENTRED    '
	#print '##################'

	def P_Roff(Roff):
		return abs((Roff/s_off**2)*np.exp(-0.5*(Roff/s_off)**2))
	
	def monopole_off(R,theta):
		argumento = lambda x: monopole(R**2+x**2-2*x*R*np.cos(theta))*P_Roff(x)
		integral1  = integrate.quad(argumento, -1.*np.inf, 0)[0]
		integral2  = integrate.quad(argumento, 0., R)[0]
		integral3  = integrate.quad(argumento, R, np.inf)[0]
		return integral1 + integral2 + integral3
	vec_moff = np.vectorize(monopole_off)

	def Delta_Sigma_off(R,theta):
		argumento = lambda x: monopole_off(x,theta)*x
		integral  = integrate.quad(argumento, 0, R)[0]
		DS_off    = (2./R**2)*integral - monopole_off(R,theta)
		return DS_off

	vec_DSoff = np.vectorize(Delta_Sigma_off)
	
	def quadrupole_off(R,theta):
		def q_off(roff):
			return quadrupole(R**2+roff**2-2*roff*R*np.cos(theta))*P_Roff(roff)
		argumento = lambda x: q_off(x)
		integral1  = integrate.quad(argumento, -1.*np.inf, 0)[0]
		integral2  = integrate.quad(argumento, 0., R)[0]
		integral3  = integrate.quad(argumento, R, np.inf)[0]
		return integral1 + integral2 + integral3	
	vec_qoff = np.vectorize(quadrupole_off)
	
	
	def psi2_off(R,theta):
		argumento = lambda x: (x**3)*monopole_off(x,theta)
		integral = integrate.quad(argumento, 0, R)[0]
		return integral*(-2./(R**2))
		
	vecpsi2_off = np.vectorize(psi2_off)	

	def quantities_centred(r):
		
		# optimize using unique r
		r,c = np.unique(r,return_counts=True)
		
		monopole_r = monopole(r)
		quadrupole_r = quadrupole(r)
		print 'computing psi2 centred'	
		psi2_r = vecpsi2(r)
		
		monopole_r   = np.repeat(monopole_r,c)
		quadrupole_r = np.repeat(quadrupole_r,c)
		psi2_r       = np.repeat(psi2_r,c)
		
		return monopole_r,quadrupole_r,psi2_r
		
	def quantities_misscentred(r):
		print 'computing misscentred profile'
		r,c = np.unique(r,return_counts=True)
		gamma_t0_off = []
		gamma_t_off = []
		gamma_x_off = []
		for R in r:
			print 'computing DS_t0_off'
			t1 = time.time()
			def DS_t_off0(theta):
				gamma_t0 = Delta_Sigma_off(R,theta)
				return gamma_t0
			argumento = lambda x: DS_t_off0(x)
			integral  = integrate.quad(argumento, 0, 2.*np.pi,points=[np.pi])[0]
			gamma_t0_off = np.append(gamma_t0_off,integral/(2.*np.pi))
			t2 = time.time()
			print (t2-t1)/60.
			
			print 'computing DS_t_off'
			t1 = time.time()
			def DS_t_off(theta):
				gamma_t0 = Delta_Sigma_off(R,theta)
				gamma_t2 = ((-6*psi2_off(R,theta)/R**2) 
				            - 2.*monopole_off(R,theta) 
				            + quadrupole_off(R,theta))
				return gamma_t0 + ellip*gamma_t2*np.cos(2.*theta)

			argumento = lambda x: DS_t_off(x)*np.cos(2.*x)
			integral  = integrate.quad(argumento, 0, 2.*np.pi,epsabs=1.e-6,epsrel=1.e-6)[0]
			gamma_t_off = np.append(gamma_t_off,integral/np.pi)
			t2 = time.time()
			print (t2-t1)/60.
			
			print 'computing DS_x_off'
			t1 = time.time()
			def DS_x_off(theta):
				gamma_x2 = ((-6*psi2_off(R,theta)/R**2) 
				            - 4.*monopole_off(R,theta))
				return ellip*gamma_x2*np.sin(2.*theta)
			argumento = lambda x: DS_x_off(x)*np.cos(2.*x)
			integral  = integrate.quad(argumento, 0, 2.*np.pi)[0]
			gamma_x_off = np.append(gamma_x_off,integral/np.pi)	
			t2 = time.time()
			print (t2-t1)/60.
			
		gamma_t0_off = np.repeat(gamma_t_off,c)
		gamma_t_off = np.repeat(gamma_t_off,c)
		gamma_x_off = np.repeat(gamma_t_off,c)
		
		return [gamma_t0_off, gamma_t_off, gamma_x_off]
		
	
	m,q,p2 = quantities_centred(r)
	output = {'monopole':m,'quadrupole':q,'psi2':p2}
	
	if misscentred:
		slicer = int(round(len(r)/ncores, 0))
		slices = ((np.arange(ncores-1)+1)*slicer).astype(int)
		pool = Pool(processes=(ncores))
		
		salida=np.array(pool.map(quantities_misscentred, np.split(r,slices))
		
		gt0_off = salida[0,0]
		gt2_off = salida[0,1]
		gx2_off = salida[0,2]
		
		for j in np.arange(1,ncores): 
			gt0_off = np.concatenate((gt0_off,salida[j,0]))
			gt2_off = np.concatenate((gt2_off,salida[j,1]))
			gx2_off = np.concatenate((gx2_off,salida[j,2]))
		
		output.update({'gt0_off':gt0_off,'gt2_off':gt2_off,'gx2_off':gx2_off})
		
	return output
