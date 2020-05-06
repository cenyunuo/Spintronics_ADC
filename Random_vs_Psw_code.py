import numpy as np
import matplotlib.pyplot as plt
from random import random
plt.rcParams.update({'font.size':15})

def comparator(In_1,In_2):
	if In_1<In_2:
		out=0
	if In_1>=In_2:
    		out=1
	
	return out

#################### constant parameters ############################
gamma=1.76e11;           # Gyromagnetic ratio [(rad)/(s.T)]
mu0=4*np.pi*1e-7 ;      # in T.m/A

q=1.6e-19;               # in Coulomb
hbar=1.054e-34;          # Reduced Planck's constant (J-s)
K_B=1.38064852e-23    #in J/K
#################### parameters related to nanomagnet ################
alpha=0.01              # Gilbert damping parameter
Ms = 1.23e6            # in A/m
d=40e-9
t=1e-9
Delta=43
#################################################################
V=(np.pi/4.0)*(d*d*t)
T=30
Hk=(2*Delta*K_B*T)/(mu0*Ms*V)
eta=0.6
Isc=((2*q*alpha)/(hbar*eta))*(mu0*(Hk)*Ms*V)
print('Isc = ' +str(Isc))
#################################################################

################### Triangular Wave details ##########################
A=2.5*Isc        # Amplitude
T=4e-9             # Time period
n_sample=3001
#################################################################
n=n_sample
t=np.linspace(0,T,n)
signal=np.zeros(n)
Total_signal=[]
Total_time=[]
########### Creating Triangular wave for a single time period ################
for i in range(len(t)):
	if t[i]<T/4.0:
		signal[i]=(4*A/T)*t[i]
	if t[i]>=T/4.0 and t[i]<(3*T/4.0):
		signal[i]=(-4*A/T)*t[i]+2*A
	if t[i]>=(3*T/4.0) and t[i]<=T:
		signal[i]=(4*A/T)*t[i]-4*A
####################################################################
cycle=2
############ Repeating the wave for the given number of cycles ##############
for rep in range(cycle):
	Total_signal=np.append(Total_signal, signal)
	Total_time=np.append(Total_time,t+rep*T)
####################################################################

Total_length=len(Total_time)
Crit_value=Isc
Crit_value=Crit_value*np.ones(Total_length) # creating critical value array


########## Probability switching calculation #################
t0=1e-9
pulse_width=T/2.0
beta=1.0
Psw=np.zeros(Total_length)
Error=np.zeros(Total_length)
for ti in range(Total_length):
	i=Total_signal[ti]/Isc
	f1=np.exp(-Delta*(1-i)**beta)
	f2=(Total_time[ti]/t0)*f1
	#f2=(pulse_width/t0)*f1
	Psw[ti]=1.0-np.exp(-f2)
	if i<=1:
		Error[ti]=abs(Psw[ti]-0.0)
	else:
		Error[ti]=abs(1.0-Psw[ti])

max_error=max(Error)
min_error=min(Error)

print('Max error = ' + str(max_error))
print('Min error = ' + str(min_error))

error_index=np.asarray(np.where(Error==max_error))
print('Error_index = ' + str(error_index[0,0]))

#index=error_index[0,0]

Psw_at_error_index = Psw[error_index]

print('Psw at max error index = ' + str(np.asarray(Psw_at_error_index)))
print('I at error index = ' + str(Total_signal[error_index]))
#exit()
################### Random switching Probability #####################
random_probability=np.zeros(Total_length)
for ti in range(Total_length):
	if Total_signal[ti]>Crit_value[ti]:
		random_probability[ti]=random()
#print(max(random_probability))
############# Comparison of Random switching vs Psw #################
compare_result=np.zeros(Total_length)
for ti in range(Total_length):
	compare_result[ti]=comparator(random_probability[ti],Psw[ti])

###################### Plot the signal #################################


Total_time=Total_time*1e9
Total_signal = Total_signal*1e6
Crit_value = Crit_value*1e6



fig = plt.figure(figsize=(16,14))

plt.subplot(2,2,1)
plt.plot(Total_time,Total_signal,linewidth = 2.0, label='Signal')
plt.plot(Total_time,Crit_value, linewidth = 2.0, label='Crit Value')
#plt.xlabel('Time(ns)')
plt.ylabel('Current ($\mu$A) ')
plt.legend(loc='upper right')
plt.grid()
'''
plt.subplot(3,1,2)
plt.plot(Total_time,Psw, linewidth=2.0)
plt.xlabel('Time(ns)')
plt.ylabel('$P_{sw}$ ')
plt.grid()
'''
plt.subplot(2,2,2)
plt.plot(Total_time,random_probability, 'b', label='$N_{rand}$')
plt.plot(Total_time,Psw, 'r', label='$P_{sw}$')
#plt.xlabel('Time(ns)')
plt.ylabel(' Probability')
plt.legend(loc='upper right')
plt.grid()

plt.subplot(2,2,3)
plt.plot(Total_time,compare_result)
plt.xlabel('Time(ns)')
plt.ylabel(' Comparison')
plt.grid()
#plt.savefig('Comparison_of_Psw_Nrand_eta_0pt6_10us.png')	

plt.subplot(2,2,4)
plt.plot(Total_time,Error,'*-')
plt.xlabel('Time(ns)')
plt.ylabel(' Error')
plt.yscale("log")
plt.grid()
plt.show()

plt.show()
