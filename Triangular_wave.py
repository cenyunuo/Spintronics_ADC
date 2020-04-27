import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':25})

A=input('Enter amplitude : ')
T=input('Enter Time period (ms) : ')
T=float(T)
f=10**3/T 
print('Frequency = ' + str(f) + ' kHz')
n=2001
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
cycle=input('Enter number of cycle of the input signal : ') 
############ Repeating the wave for the given number of cycles ##############
for rep in range(cycle):
	Total_signal=np.append(Total_signal, signal)
	Total_time=np.append(Total_time,t+rep*T)
####################################################################
Total_length=len(Total_time)
Crit_value=input('Enter critical value : ')
Crit_value=Crit_value*np.ones(Total_length) # creating critical value array
###################### Plot the signal #################################
fig = plt.figure(figsize=(14,9))
plt.plot(Total_time,Total_signal,linewidth = 3.0, label='Signal')
plt.plot(Total_time,Crit_value, linewidth = 3.0, label='Crit Value')
plt.xlabel('Time(ms)')
plt.ylabel('Signal ')
plt.legend(loc='upper right')
plt.grid()
plt.show()
		
