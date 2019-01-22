from functions import *
import os
import matplotlib.pyplot as plt
import sys

#Define variables
searchfor="'total energy'"
searchfor2="'Time for block'"
imin=1
if len(sys.argv) > 1:
   filename=sys.argv[1]
else:
   filename='../he4/he4lin.dmc10000.out'
N=100 #number of bunching intervals
fs=17
uncertainty=[]
version=1 #0=old, 1=new

os.system("grep "+searchfor+" "+filename+" > temp.dat")
os.system("grep "+searchfor2+" "+filename+" > temp2.dat")
 
#Get data and determine uncertainties for blocks
rawmean=[]
rawerr=[]
f = open("temp.dat")
for line in f.readlines():
   rawmean.append(float(line[47:64]))
   rawerr.append(float(line[71:-4]))
f.close()
mdata=rawmean[1:len(rawmean)]
edata=rawerr[1:len(rawerr)]
os.system("rm temp.dat")

blocktime=[]
f = open("temp2.dat")
for line in f.readlines():
   blocktime.append(float(line[40:50]))
f.close()
os.system("rm temp2.dat")

#calculate and print average block time
print 'Average Block Time = '+str(sum(blocktime)/len(blocktime))+' sec'

mean=average(mdata)
print 'E = '+str(mean)+' MeV'

if version==0:
   N=5
   interval=[1]
   length=len(mdata)
   for n in range(1,N+1):
      if n>1:
         interval.append(interval[n-2]*2)
      if length%2 == 0:
         odd=False
      else:
         odd=True
      tempmean=[]
      uncertainty.append(sig(mdata))
      length=length/2
      for i in range(0,length/2-1):
         if (i==length/2-1 and odd):
            tempmean.append(1.0/3.0*(mdata[2*i]+mdata[2*i+1]+mdata[2*i+2]))
         else:
            tempmean.append(0.5*(mdata[2*i]+mdata[2*i+1]))
      mdata=tempmean
elif version==1:
   interval=[]
   length=len(mdata)
   for nblk in range(1,N+1): #nblk is the block size (number being averaged together)
      tempmean=[]
      nblks=length/nblk # nblks is the number of blocks for each block size
      interval.append(nblk)
      idx=imin-1
      for i in range(0,nblks):
         temp=0.0
         for j in range(0,nblk):
            temp=temp+mdata[idx]
            idx=idx+1
         tempmean.append(temp/nblk)
      uncertainty.append(sig(tempmean))
else:
   print 'Error: version must be 0 (old) or 1 (new)'

#print interval
#print uncertainty
f1 = plt.figure()
ax1 = f1.add_subplot(111)
f1.suptitle('Binding energy', fontsize=16)
plt.xlabel('ave block')
plt.ylabel('Energy (MeV)')
ax1.errorbar(range(1,len(rawmean)+1),rawmean,yerr=rawerr,fmt='ro')

f2 = plt.figure()
ax2 = f2.add_subplot(111)
f2.suptitle('Error as a function of block size', fontsize=16)
plt.xlabel('block size')
plt.ylabel('error')
plt.xlim([min(interval),max(interval)])
ax2.plot(interval, uncertainty,'bo-')
plt.show()
