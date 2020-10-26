import numpy as np
from random import gauss
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d 
plt.rcParams.update({'font.size':15})
#################### Real unit  LLGS ###########################
def LLGS(m,H,mp,Beta):
	global alpha,mu0,gamma

	precision=-((gamma*mu0)/(1+alpha**2))*np.cross(m,H)
	damping=-((gamma*mu0*alpha)/(1+alpha**2))*np.cross(m,np.cross(m,H))
	slonczewski=-(Beta/(1+alpha**2))*np.cross(m,np.cross(m,mp))
	field_like=((alpha*Beta)/(1+alpha**2))*np.cross(m,mp)

	dmdt=precision+damping+slonczewski+field_like

	return dmdt
#################################################################

##################### Vector magnitude #############################
def mag(M):
	magnitude = np.sqrt((M[0])**2+(M[1])**2+(M[2])**2)
	return magnitude
#################################################################

##################### sec hyparabolic function #######################
def sech(x):
	y=2.0/(np.exp(x)+np.exp(-x))
	return(y)
################### LLGS calculation result #########################
def calculate_LLGS(m,Hk,Hd,mp,Beta):
	global n, alpha, K_B, gamma, Ms, V, delta_t, mu0, T
	h_step = delta_t
	ti=0
	m_tilda=np.zeros((n,3))
	H_uni=np.zeros(3)
	H_demag=np.zeros(3)
	
	while ti<(n-1):
		#print(ti)
		mx=m[ti,0]
		my=m[ti,1]
		mz=m[ti,2]
		############## Uniaxial Field ######################
		H_uni=np.array([0, 0, Hk*mz])
		############## Demag Field #######################
		H_demag=np.array([Hd[0]*mx,Hd[1]*my,Hd[2]*mz])
		############### Thermal section ###################
		G01=np.zeros(3)
		G01[0]=gauss(0,1)
		G01[1]=gauss(0,1)
		G01[2]=gauss(0,1)
		alpha_const=alpha/(1+alpha**2)
	
		T_K=0.0*T     # in Kelvin
		const=(2*K_B*T_K)/(gamma*Ms*V*delta_t)  
		Therm_const=np.sqrt(alpha_const*const)/mu0    # unit in A/m
		H_therm=Therm_const*G01

		H_eff=H_uni+H_demag+H_therm

		m_tilda[ti+1,:]=m[ti,:]+h_step*LLGS(m[ti,:],H_eff,mp,Beta)
		m[ti+1,:]=m[ti,:]+(h_step/2.0)*(LLGS(m[ti,:],H_eff,mp,Beta)+LLGS(m_tilda[ti+1,:],H_eff,mp,Beta))

		ti=ti+1

	return m

######################## SHM Class ###############################
class SHM:	
	def __init__(self, W, L, t):
	     self.W=W
	     self.L=L
	     self.t=t

#################### constant parameters ############################
gamma=1.76e11;           # Gyromagnetic ratio [(rad)/(s.T)]
mu0=4*np.pi*1e-7 ;      # in T.m/A

q=1.6e-19;               # in Coulomb
hbar=1.054e-34;          # Reduced Planck's constant (J-s)
K_B=1.38064852e-23    #in J/K
T=300                             # in Kelvin
#################### parameters related to nanomagnet ################
alpha=0.007              # Gilbert damping parameter
Ms=850*1e3            # in A/m
Delta=42

Length=120e-9
Width=40e-9
t_FL=1.5e-9

A_MTJ=(np.pi/4.0)*Length*Width 
Area=A_MTJ

lambda_sf=1.4e-9            # spin flip length

V=A_MTJ*t_FL;             # Volume [m^3]
E_B=Delta*K_B*T

Ku2=E_B/V              # in J/m^3
################# Anisotropy field ###################################
Hk=2*Ku2/(mu0*Ms)             # Uniaxial field [in A/m]

################## Demagnetization field related #######################

Nxx=0.066
Nyy=0.911
Nzz=0.022

hdx=-Nxx*Ms
hdy=-Nyy*Ms
hdz=-Nzz*Ms
Hd=np.array([hdx, hdy, hdz])
######################## SHM parameter values #######################
##### SHM(W,L,t)
SHM_1 = SHM(120e-9,80e-9,2.1e-9)
#SHM_2 = SHM(150e-9,100e-9,2.8e-9)

SHM_2 = SHM(200e-9,100e-9,2.5*2.8e-9)


A_SHM_1=SHM_1.W*SHM_1.t
A_SHM_2=SHM_2.W*SHM_2.t

theta_SHE=0.3

###################################################################

P_she_1=(A_MTJ/A_SHM_1)*theta_SHE*(1.0-sech(SHM_1.t/lambda_sf))
P_she_2=(A_MTJ/A_SHM_2)*theta_SHE*(1.0-sech(SHM_2.t/lambda_sf))

Ic = float(input('Enter charge current(in uA) = '))
Ic=Ic*1e-6

###################### time related portion ###########################
stop_time=2e-9             # in s
n=50001
t=np.linspace(0,stop_time,n) # in ns
h_step=t[1]-t[0]
delta_t=h_step
##################### Initial Magnetization ##########################
m_1=np.zeros((n,3))            
m_2=np.zeros((n,3))

mz0=-0.9999
mx0=np.sqrt(1-mz0**2)
my0=0

m_1[0,:]=[mx0,my0,mz0]
m_2[0,:]=[mx0,my0,mz0]

#################### Spin current section ############################
mp=[0,0.0,1.0]           # dimensionless spin current vector

Is_1 = P_she_1*Ic           # spin current entering into MTJ 1
Is_2 = P_she_2*Ic           # spin current entering into MTJ 2
print('Is_1 = %e' %Is_1)
print('Is_2 = %e' %Is_2)

#exit()

J_MTJ_1=Is_1/A_MTJ
J_MTJ_2=Is_2/A_MTJ

Beta_1= (gamma*hbar*J_MTJ_1)/(2*q*Ms*t_FL)
Beta_2= (gamma*hbar*J_MTJ_2)/(2*q*Ms*t_FL)

print('MTJ 1 is being simulated...')
M_1 = calculate_LLGS(m_1,Hk,Hd,mp,Beta_1)
print('MTJ 1 simulation is done')

print('MTJ 2 is being simulated...')
M_2 = calculate_LLGS(m_2,Hk,Hd,mp,Beta_2)
print('MTJ 2 simulation is done')

t=t*1e9 # time converted to ns scale

fig = plt.figure(figsize=(14,7))

plt.subplot(1,2,1)
plt.plot(t,M_1[:,2])
plt.grid()
plt.xlabel('Time(ns)')
plt.ylabel(r"$m_z$")
plt.title(r"$MTJ_L$")
plt.ylim([-1.2,1.2])

plt.subplot(1,2,2)
plt.plot(t,M_2[:,2])
plt.grid()
plt.xlabel('Time(ns)')
plt.ylabel(r"$m_z$")
plt.title(r"$MTJ_H$")
plt.ylim([-1.2,1.2])
plt.show()
