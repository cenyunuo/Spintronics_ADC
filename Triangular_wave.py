import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':15})

A=input('Enter amplitude : ')
T=input('Enter Time period : ')
T=float(T)
n=2001
t=np.linspace(0,T,n)
signal=np.zeros(n)
Total_signal=[]
Total_time=[]
for i in range(len(t)):
	if t[i]<T/4.0:
		signal[i]=(4*A/T)*t[i]
	if t[i]>=T/4.0 and t[i]<(3*T/4.0):
		signal[i]=(-4*A/T)*t[i]+2*A
	if t[i]>=(3*T/4.0) and t[i]<=T:
		signal[i]=(4*A/T)*t[i]-4*A
cycle=input('Enter number of cycle of the input signal : ')
for rep in range(cycle):
	Total_signal=np.append(Total_signal, signal)
	Total_time=np.append(Total_time,t+rep*T)

Total_length=len(Total_time)
Crit_value=input('Enter critical value : ')
Crit_value=Crit_value*np.ones(Total_length)
'''
plt.subplot(2,1,1)
plt.plot(t,signal)
plt.grid()
'''
#plt.subplot(2,1,2)
plt.plot(Total_time,Total_signal)
plt.plot(Total_time,Crit_value)
plt.grid()
plt.show()
		
