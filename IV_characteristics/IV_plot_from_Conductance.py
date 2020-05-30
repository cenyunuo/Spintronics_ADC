import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':20})
import xlsxwriter 
def Gcal(c1,c2,c3,c4,t,V):
	G=np.exp(c1*t+c2)*V**2+np.exp(c3*t+c4)
	return G
################# Constant values #########################
c1_p=-5.90497
c2_p=21.5434
c3_p=-7.46919
c4_p=25.0243

c1_ap=-6.7524
c2_ap=23.2848
c3_ap=-7.56891
c4_ap=24.144
################ Parameter values #########################
Tox=1.15 # should be in nm
t_MgO=Tox
W=200e-9
L=100e-9
Area=(W*L*1e4)
########################################################
nV=51
V=np.linspace(0,3,nV)
Ip=np.zeros(nV)
Iap=np.zeros(nV)
Rp=np.zeros(nV)
Rap=np.zeros(nV)
########################################################
workbook = xlsxwriter.Workbook('I-V_result.xlsx') 
worksheet = workbook.add_worksheet() 
worksheet.write('A1', 'Voltage') 
worksheet.write('B1', 'Ip') 
worksheet.write('C1', 'Rp') 
worksheet.write('D1', 'Iap') 
worksheet.write('E1', 'Rap') 

############### I-V characteristics for parallel ###############
for i in range(nV):
	worksheet.write(i+1,0,V[i])
	G=Gcal(c1_p,c2_p,c3_p,c4_p,t_MgO,V[i]) # conductance per unit area
	G=G*Area
	Rp[i]=1/G
	Ip[i]=V[i]/Rp[i]
	worksheet.write(i+1,1,Ip[i])
	worksheet.write(i+1,2,Rp[i])
############### I-V characteristics for Antiparallel ###############
for i in range(nV):
	G=Gcal(c1_ap,c2_ap,c3_ap,c4_ap,t_MgO,V[i]) # conductance per unit area
	G=G*Area
	Rap[i]=1/G
	Iap[i]=V[i]/Rap[i]
	worksheet.write(i+1,3,Iap[i])
	worksheet.write(i+1,4,Rap[i])
###########################################################
workbook.close() 

fig = plt.figure(figsize=(12,9))

plt.subplot(2,1,1)
plt.plot(V,Ip)
plt.grid()
plt.xlabel('V')
plt.ylabel('$I_p$')
#plt.title('Parallel')

plt.subplot(2,1,2)
plt.plot(V,Iap)
plt.grid()
plt.xlabel('V')
plt.ylabel('$I_{ap}$')
#plt.title('Antiparallel')

plt.show()
