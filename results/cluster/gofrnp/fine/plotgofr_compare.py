from numpy import *
#import numpy as np
import matplotlib.pyplot as plt

#set up input parameters
calc=2 #1=14n2p 2=2n2p
corrtype=1 #1=lin, 2=ip
sys=5 #np0=1, np1=3, pp=5, nn=7
#filenames=('temp.dmc','gofrnp_0.00025_lin_alpha.dmc','gofrnp_0.00025_ip_alpha.dmc','gofrnp_0.0005_lin.dmc','gofrnp_0.0005_ip.dmc')
#titles=('Continuum Alpha','0.00025 lin Alpha','0.00025 ip Alpha','0.00025 lin Cluster','0.00025 ip Cluster')
#pltstyle=('ko','ro','mo','yo','go','bo')
filenames=('../gofrnp_continuum_ip.dmc','gofrnp_0.00025_lin_fine.dmc','gofrnp_0.00025_ip_fine.dmc')
titles=('Alpha Continuum','Linear - 0.00025','IP - 0.00025')
#filenames=('../gofrnp_continuum_ip.dmc','gofrnp_0.00025_lin_alpha_fine.dmc','gofrnp_0.00025_ip_alpha_fine.dmc')
#titles=('Alpha Continuum','Linear 2n2p - 0.00025','IP 2n2p - 0.00025')
pltstyle=('ko','rs','b^')

#build other parameters
lenfiles=len(filenames)
if(corrtype==1):
   corr='lin'
elif(corrtype==2):
   corr='ip'
else:
   print ''
   print 'Invalid choice for corrtype! Valid choices are 1 and 2. Exiting!'
   print ''
   exit()
#if(not(sys==1 or sys==3 or sys==5 or sys==7)):
#   print 'Invalid choice for system! Valid choices are 1, 3, 5, and 7. Exiting!'
#   exit()
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
if(calc==1):
   calcname=''
elif(calc==2):
   calcname='_alpha'
else:
   print ''
   print 'Invalid choice for calculation! Valid choices are 1 and 2. Exiting!'
   print ''
   exit()

data=loadtxt(filenames[0])
for n in range(0,lenfiles):
   data=array(loadtxt(filenames[n]))
#   print data[:,0]
#   data[:,0] = (data[:,0])[0::2]
#   data[:,0] = delete(data[:,0], arange(0, data[:,0].size, 2))
#   print ''
#   print data[:,0]
#   stop

   #normalize data to have unit integral
   norm=0.0
   length=len(data)
   for i in range(0,length):
      norm=norm+data[i,sys]*(data[1,0]-data[0,0])
   for i in range(0,length):
      data[i,sys]=data[i,sys]/norm
      data[i,sys+1]=data[i,sys+1]/norm

   #Plot the data
   f1=plt.figure(1)
   plt.rcParams.update({'font.size': 15})
   plt.xlim([0,10])
   plt.xlabel('r(fm)')
   plt.ylabel(r'g$_{'+sysname+'}(r)$')
   plt.axhline(linewidth=2, color='k')
   if(n==0):
      plt.errorbar((data[:,0])[0::3],(data[:,sys])[0::3],yerr=(data[:,sys+1])[0::3],fmt=pltstyle[n],label=titles[n])
   else:
      plt.errorbar(data[:,0],data[:,sys],yerr=data[:,sys+1],fmt=pltstyle[n],label=titles[n])
   plt.legend(loc='upper right',numpoints=1)

#   f2=plt.figure(2)
#   plt.xlim([4,10])
#   plt.ylim([0,0.05])
#   plt.xlabel('r(fm)')
#   plt.ylabel(r'g$_{'+sysname+'}(r)$')
#   plt.errorbar(data[n][:][0],data[n][:][sys],yerr=data[n][:][sys+1],fmt=pltstyle[n],label=titles[n])
#   plt.legend(loc='upper right',numpoints=1)

f1.savefig("plot1.pdf", bbox_inches='tight')
#f2.savefig("plot2.pdf", bbox_inches='tight')
plt.show()
