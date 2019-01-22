from functions import *
import os
import matplotlib.pyplot as plt
import sys
import numpy as np

#set up parameters
imin=1
tot=100 #total number of files

#get data
data=[]
tdata=[]
for i in range(imin,tot+1):
   new=[]
   filename="he4lin.unc10000_"+str(i)+".out"
   os.system("grep 'total energy' "+filename+" > temp.dat")
   f = open("temp.dat")
   for line in f.readlines():
      new.append(float(line[47:64]))
      if i==imin:
         tdata.append(float(line[31:45]))
   f.close()
   os.system("rm temp.dat")
   data.append(new)

#average data
avedata=[]
errdata=[]
length=len(data[:][0])
print data
if tdata==[]:
   tdata=range(1,length+1)
for i in range(0,length):
   total=[]
   for j in range(imin-1,tot):
      total.append(data[j][i])
   avedata.append(average(total))
   errdata.append(sig(total))
#   avedata.append(np.mean(total))
#   errdata.append(np.std(total))
print avedata
print errdata
plt.errorbar(tdata,avedata,yerr=errdata,fmt='o',label="linear")
plt.legend()
plt.title("Unconstrained calculation")
plt.xlabel(r'$\tau$')
plt.ylabel(r'E')
plt.ylim([-30,-24])
plt.show()
