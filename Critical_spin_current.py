import numpy as np
########## Constant parameters ########
q=1.60217662e-19    # in coulomb
hbar=1.054e-34;          # Reduced Planck's constant (J-s)
alpha=0.01
K_B = 1.38064852e-23    #in J/K
mu0 = 4*np.pi * 1.0e-7      #in N/A^2
########## MTJ parameter values ########
Ms=450e3 # in A/m
Delta=45
T=300 

Major_length=input('Enter Major axis length (nm)= ')
a=(float(Major_length)/2.0)*1e-9
Minor_length=input('Enter Minor axis length (nm)= ')
b=(float(Minor_length)/2.0)*1e-9
Thickness=1e-9
Area=np.pi*a*b
Volume=Area*Thickness

Ku=(Delta*K_B*T)/Volume
Hk=(2*Ku)/(mu0*Ms)
Hd=Ms

Heff=abs(Hk-Hd)
Ic=(2*q*alpha/hbar)*(mu0*Heff*Ms*Volume)

print('Critical current  = ' + str(Ic))
