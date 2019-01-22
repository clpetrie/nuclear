from functions import *
import os
import matplotlib.pyplot as plt
import sys
import numpy as np

#set up parameters
imin=1
tot=100 #total number of files
file1="lindata_fine/he4lin.unc10000_"
file2="ipdata_fine/he4ip.unc10000_"

#get data
data1=[]
tdata1=[]
for i in range(imin,tot+1):
   new=[]
   filename=file1+str(i)+".out"
   os.system("grep 'total energy' "+filename+" > temp.dat")
   f = open("temp.dat")
   for line in f.readlines():
      new.append(float(line[47:64]))
      if i==imin:
         tdata1.append(float(line[31:45]))
   f.close()
   os.system("rm temp.dat")
   data1.append(new)
data2=[]
tdata2=[]
for i in range(imin,tot+1):
   new=[]
   filename=file2+str(i)+".out"
   os.system("grep 'total energy' "+filename+" > temp.dat")
   f = open("temp.dat")
   for line in f.readlines():
      new.append(float(line[47:64]))
      if i==imin:
         tdata2.append(float(line[31:45]))
   f.close()
   os.system("rm temp.dat")
   data2.append(new)


#average data
avedata1=[]
errdata1=[]
length=len(data1[:][0])
if tdata1==[]:
   tdata1=range(1,length+1)
for i in range(0,length):
   total=[]
   for j in range(imin-1,tot):
      total.append(data1[j][i])
   avedata1.append(average(total))
   errdata1.append(sig(total))
avedata2=[]
errdata2=[]
length=len(data2[:][0])
if tdata2==[]:
   tdata2=range(1,length+1)
for i in range(0,length):
   total=[]
   for j in range(imin-1,tot):
      total.append(data2[j][i])
   avedata2.append(average(total))
   errdata2.append(sig(total))

plt.errorbar(tdata1,avedata1,yerr=errdata1,fmt='bo',label="linear")
plt.errorbar(tdata2,avedata2,yerr=errdata2,fmt='ro',label="independent pair")
plt.legend()
plt.title("Unconstrained calculation")
plt.xlabel(r'$\tau$')
plt.ylabel(r'E')
if file1=="lindata/he4lin.unc10000_":
   plt.xlim([0,0.01])
plt.ylim([-30,-24])
#plt.ylim([-28,-26])
plt.show()

