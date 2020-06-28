import numpy as np
from random import gauss
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d 
plt.rcParams.update({'font.size':18})

##################### Effective Magnetic field #####################
def H_eff(m, H_ext):
	global Hk, hdx, hdy, hdz
	mx=m[0]
	my=m[1]
	mz=m[2]
	H_uni=np.array([0, 0, Hk*mz])
	H_demag=np.array([hdx*mx, hdy*my, hdz*mz])
	H_effective=H_uni+H_demag+H_ext

	return H_effective

#################### Real unit  LLGS ###########################
def LLGS(m,H,Zeta_SHE,mp):
	global alpha,mu0,gamma

	precision=-((gamma*mu0)/(1+alpha**2))*np.cross(m,H)
	damping=-((gamma*mu0*alpha)/(1+alpha**2))*np.cross(m,np.cross(m,H))
	slonczewski=-(Zeta_SHE/(1+alpha**2))*np.cross(m,np.cross(m,mp))
	field_like=((alpha*Zeta_SHE)/(1+alpha**2))*np.cross(m,mp)

	dmdt=precision+damping+slonczewski#+field_like

	return dmdt
##############################################################
##################### Vector magnitude #############################
def mag(M):
	magnitude = np.sqrt((M[0])**2+(M[1])**2+(M[2])**2)
	return magnitude
#################################################################

#################### constant parameters ############################
gamma=1.76e11;           # Gyromagnetic ratio [(rad)/(s.T)]
mu0=4*np.pi*1e-7 ;      # in T.m/A

q=1.6e-19;               # in Coulomb
hbar=1.054e-34;          # Reduced Planck's constant (J-s)
K_B=1.38064852e-23    #in J/K
#################### parameters related to nanomagnet ################
alpha=0.01              # Gilbert damping parameter
Ms=450e3               # in A/m

t_FL=1.0e-9
Length=60e-9
Width=60e-9
Area=(Length*Width) * (np.pi/4.0)
A_MTJ=Area
magVolume = t_FL *  Area  # in m^3
V=magVolume
################ Anisotropy field related ################
delta=45
Temperature = 300              # in K

Eb=delta*K_B*Temperature
Ku2=Eb/magVolume+(mu0*Ms**2*0.5)

#print(Ku2)
#Ku2=2.245e5       # in J/m^3
Hk=(2*Ku2)/(mu0*Ms)
print('Anisotropy field magnitude = ' + str(mu0*Hk) + ' Tesla')

H_uni=np.zeros(3)
############ Demagnetization field related ###############

Nxx=0*0.0123
Nyy=0*0.0123
Nzz=1.0-(Nxx+Nyy)

hdx=-Nxx*Ms
hdy=-Nyy*Ms
hdz=-Nzz*Ms
print('Demag magnitude along -ve z = ' + str(mu0*Ms) + ' Tesla')
H_demag=np.zeros(3)
#exit()
############### External Mag field #####################
H_ext_mag=0.0*38e-3 # in Tesla along y gives the small oscillation

H_ext_mag=H_ext_mag/mu0
H_ext=np.array([H_ext_mag,0,0])

Ic=20e-6
J=Ic/A_MTJ
Zeta_SHE=(gamma*hbar*J)/(2*q*t_FL*Ms)
mp=np.array([0,0,-1])
############ time related portion #######################
stop_time=20.0e-9             # in ns
n=20001
t=np.linspace(0,stop_time,n) # in ns
h_step=t[2]-t[1]
delta_t=h_step
print(delta_t)
#exit()
############# Initial Magnetization #####################
m=np.zeros((n,3))
mz0=0.99
mx0=np.sqrt(1-mz0**2)
my0=0
m[0,:]=[mx0,my0,mz0]
print(m[0,:])

ti=0
m_tilda=np.zeros((n,3))
m_mag=np.zeros(n)
m_mag[ti]=mag(m[ti,:])


while ti<(n-1):
	print('-------------------------------------------------------------------------------------------')
	
	k1=LLGS(m[ti,:],H_eff(m[ti,:],H_ext),Zeta_SHE,mp)
	#print(h_step*k1/2.0)
	m_k1=(m[ti,:]+h_step*k1/2.0)
	k2=LLGS(m_k1,H_eff(m_k1,H_ext),Zeta_SHE,mp)  
	#print(h_step*k2/2.0)
	m_k2= (m[ti,:]+h_step*k2/2.0)     
	k3=LLGS(m_k2,H_eff(m_k2,H_ext),Zeta_SHE,mp)
	#print(h_step*k3)
	m_k3=(m[ti,:]+h_step*k3)
	k4=LLGS(m_k3,H_eff(m_k3,H_ext),Zeta_SHE,mp)
	
	m[ti+1,:]=m[ti,:]+(h_step/6.0)*(k1+2*k2+2*k3+k4)

	print(ti)
	print(m[ti+1,2])
	m_mag[ti+1]=mag(m[ti+1,:])
	ti=ti+1

############## Plot the result ########################
fig = plt.figure(figsize=(15,12))

plt.subplot(2,2,1)
plt.plot(t*1e9,m[:,0],linewidth = 2.5)
plt.grid()
plt.xlabel('Time(ns)')
plt.ylabel(r"$m_x$")
plt.ylim([-1.2,1.2])
#plt.show()
	
plt.subplot(2,2,2)
plt.plot(t*1e9,m[:,1],linewidth = 2.5)
plt.grid()
plt.xlabel('Time(ns)')
plt.ylabel(r"$m_y$")
plt.ylim([-1.2,1.2])
#plt.show()

plt.subplot(2,2,3)
plt.plot(t*1e9,m[:,2],linewidth = 2.5)
plt.grid()
plt.xlabel('Time(ns)')
plt.ylabel(r"$m_z$")
plt.ylim([-1.2,1.2])


plt.subplot(2,2,4)
plt.plot(t*1e9,m_mag,linewidth = 2.5)
plt.grid()
plt.xlabel('Time(ns)')
plt.ylabel(r"$|m|$")
plt.ylim([-1.2,1.2])

plt.show()

