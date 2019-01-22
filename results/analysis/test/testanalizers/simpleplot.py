import os
import matplotlib.pyplot as plt

#Define variables
filename=['unconstrained/he4lin_vmc10000.dat','unconstrained/he4lin_dmc10000.dat','unconstrained/he4lin_unc10000.dat']
#filename=['unconstrained/he4lin_unc10000.dat']
fs=17

#calculations for no ip and print those
for s in range(0,len(filename)):
#   tdata=[]
   mdata=[]
   edata=[]
   f = open(filename[s])
   for line in f.readlines():
#      tdata.append(float(line[31:45]))
      mdata.append(float(line[47:64]))
      edata.append(float(line[71:-4]))
   f.close()

   mean=0
   err=0
   count=0
   length=len(mdata)
   for i in range(0,length):
      count=count+1
      mean=mean+mdata[i]
      err=err+edata[i]**2
   mean=mean/(count)
#   err=(err/(count))**(0.5)
   err=(err)**(0.5)
   os.system('echo E'+str(s)+' = '+str(mean)+' +/- '+str(err)+' MeV')

#   plt.plot(tdata,mdata,label='E'+str(s)+' = '+str(round(mean,2))+'+/-'+str(round(err,2))+' MeV')
   plt.plot(mdata,'-o',label='E'+str(s)+' = '+str(round(mean,2))+'+/-'+str(round(err,2))+' MeV')
#   plt.xlabel(r'$\tau$ (MeV${}^{-1}$)',fontsize=fs)
   plt.xlabel(r'Averaging Block',fontsize=fs)
   plt.ylabel('E (MeV)',fontsize=fs)
plt.legend(loc='upper center')
plt.show()
