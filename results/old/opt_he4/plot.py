import os
import platform
import sys
import csv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

#set up parameters
num=40 #number of parameters
filename='he4ip_optimization.out'

if len(sys.argv) > 1:
   filename=str(sys.argv[1])
#read in parameter data
os.system('grep "new parameters" '+filename+' > temp.txt')
params=np.loadtxt('temp.txt',usecols=(range(3,num+3)))
os.system('rm temp.txt')
#read in energy data
os.system('grep "total energy" '+filename+' > temp.txt')
energy=np.loadtxt('temp.txt',usecols=range(3,4))
eerr=np.loadtxt('temp.txt',usecols=(range(5,6)))
os.system('rm temp.txt')

#do plot of energy then parameters
f1=plt.figure()
x=range(1,len(energy)+1)
plt.close('all')
plt.errorbar(x,energy,yerr=eerr,fmt='-o')

#plt.figure()
f2, axarr = plt.subplots(14, figsize=(8,15), sharex=True)
axarr[0].plot(params[:,0])
axarr[1].plot(params[:,1])
axarr[2].plot(params[:,2])
axarr[3].plot(params[:,3])
axarr[4].plot(params[:,4])
axarr[5].plot(params[:,5])
axarr[6].plot(params[:,6])
axarr[7].plot(params[:,7])
axarr[8].plot(params[:,8])
axarr[9].plot(params[:,11])
axarr[10].plot(params[:,12])
axarr[11].plot(params[:,13])
axarr[12].plot(params[:,14])
axarr[13].plot(params[:,15])
axarr[0].set_ylabel('test')
axarr[1].set_ylabel('test')
axarr[2].set_ylabel('test')
axarr[3].set_ylabel('test')
axarr[4].set_ylabel('test')
axarr[5].set_ylabel('test')
axarr[6].set_ylabel('test')
axarr[7].set_ylabel('test')
axarr[8].set_ylabel('test')
axarr[9].set_ylabel('test')
axarr[10].set_ylabel('test')
axarr[11].set_ylabel('test')
axarr[12].set_ylabel('test')
axarr[13].set_ylabel('test')
axarr[13].set_xlabel('iteration')
plt.show()
#pp = PdfPages('plot.pdf')
#pp.savefig(f1,f2)
#pp.close()
