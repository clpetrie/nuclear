import numpy as np
import matplotlib.pyplot as plt

#set up input parameters
sys=5 #np0=1, np1=3, pp=5, nn=7
density=('0.0005','0.01','0.0005','0.01')
pltstyle=('ko','bo','go','ro')
density=('0.0005','0.0005','0.01','0.01')
pltstyle=('ko','co','bo','ro')

#build other parameters
lenfiles=len(density)
corr=('lin','ip','lin','ip')
if(sys==1):
   sysname='np0'
elif(sys==3):
   sysname='np1'
elif(sys==5):
   sysname='pp'
elif(sys==7):
   sysname='nn'
else:
   print ''
   print 'Invalid choice for system! Valid choices are 1, 3, 5, and 7. Exiting!'
   print ''
   exit()
   

data=[None]*lenfiles
for n in range(0,lenfiles):
   #build data and remove unwanted pieces
   f=open('gofrnp_'+density[n]+'_'+corr[n]+'.dmc')
   data[n]=[]
   for line in f.readlines():
      data[n].append(line.split())
   f.close()
   for i in range(0,3): #delete header information
      del data[n][0]
   length=len(data[n])-2
   for i in range(length,length+2): #delete empty two rows at the end
      del data[n][length]
   for i in range(0,length):
      for j in range(0,9): #9 is the number of systems + r + errors for each system
         data[n][i][j]=float(data[n][i][j])
   #normalize data to have unit integral
   norm=0.0
   for i in range(0,length):
      norm=norm+data[n][i][sys]*(data[n][1][0]-data[n][0][0])
   for i in range(0,length):
      data[n][i][sys]=data[n][i][sys]/norm
      data[n][i][sys+1]=data[n][i][sys+1]/norm
   data[n]=zip(*data[n])

f1=plt.figure(1)
plt.rcParams.update({'font.size': 15})
plt.xlim([0,10])
plt.xlabel('r(fm)')
plt.ylabel(r'g$_{'+sysname+'}(r)$')
plt.axhline(linewidth=2, color='k')
for n in range(0,lenfiles):
   plt.errorbar(data[n][:][0],data[n][:][sys],yerr=data[n][:][sys+1],fmt=pltstyle[n],label=r'$\rho$ = '+density[n]+r' fm$^{-3}$, '+corr[n])
plt.legend(loc='upper right',numpoints=1)
f1.savefig("plot1.pdf", bbox_inches='tight')

f2=plt.figure(2)
plt.xlim([4,10])
plt.ylim([0,0.05])
plt.xlabel('r(fm)')
plt.ylabel(r'g$_{'+sysname+'}(r)$')
for n in range(0,lenfiles):
   plt.errorbar(data[n][:][0],data[n][:][sys],yerr=data[n][:][sys+1],fmt=pltstyle[n],label=r'$\rho$ = '+density[n]+r' fm$^{-3}$, '+corr[n])
plt.legend(loc='upper right',numpoints=1)
f2.savefig("plot2.pdf", bbox_inches='tight')
plt.show()
