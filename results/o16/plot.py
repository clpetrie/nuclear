import os
import matplotlib.pyplot as plt

#Define variables
order=['Elin','Eip','Equad']
filename=['o16lin.dmc1000.dat','o16ip.dmc1000.dat','o16quad.dmc1000.dat']
col=['r','b','g']
style=['o-','s-','^-']
lab=[r'$E_{lin}$',r'$E_{ip}$',r'$E_{quad}$']
minval=[4,15,15]
fs=18 #axes font size
don=[0,1,2]
smax=len(don)



#calculations for no ip and print those
for s in range(0,smax):
   n=don[s]
   tdata=[]
   mdata=[]
   edata=[]
   f = open(filename[n])
   for line in f.readlines():
      tdata.append(float(line[31:45]))
      mdata.append(float(line[47:64]))
      edata.append(float(line[71:-4]))
   f.close()

   mean=0
   err=0
   count=0
   length=len(mdata)
   for i in range(minval[s]-1,length):
      count=count+1
      mean=mean+mdata[i]
      err=err+edata[i]**2
   mean=mean/(count)
   err=(err/(count))**(0.5)
   os.system('echo '+order[n]+' = '+str(mean)+' +/- '+str(err)+' MeV')

   plt.plot(tdata,mdata,col[n]+style[n],label=lab[n]+' = '+str(round(mean,2))+'+/-'+str(round(err,2))+' MeV')
   plt.axvline(x=tdata[minval[n]],color=col[n])
   plt.xlabel(r'$\tau$ (fm${}^{-1}$)',fontsize=fs)
   plt.ylabel('E (MeV)',fontsize=fs)
plt.legend(loc='upper center')
plt.show()
