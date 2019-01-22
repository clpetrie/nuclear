import os
import matplotlib.pyplot as plt

#Define variables
filename="../he4/he4lin_dmc10000.dat"
N=10 #max number of divisions
fs=17
uncertainty=[]

#Define some functions
def sig(indata):
   sigma=0
   length=len(indata)
   for i in range(0,length):
      sigma=sigma+indata[i]**2
   sigma=(sigma/length)**(0.5)
   return sigma

def average(indata):
   ave=0
   length=len(indata)
   for i in range(0,length):
      ave=ave+indata[i]
   ave=ave/length
   return ave

def sd(indata):
   stdev=0
   length=len(indata)
   mn=average(indata)
   for i in range(0,length):
      stdev=stdev+(indata[i]-mn)**2
   stdev=(stdev/(length-1))**(0.5)
   return stdev
 
#Get data and determine uncertainties for blocks
mdata=[]
edata=[]
f = open(filename)
for line in f.readlines():
   mdata.append(float(line[47:64]))
   edata.append(float(line[71:-4]))
f.close()

mean=average(mdata)
print 'E = '+str(mean)+' MeV'
for n in range(1,N+1):
   temparr=[]
   length=len(mdata)/n
   Npone=len(mdata)-n*length #how many have an extra (plus one)
   start=0
   end=length
   for ni in range(1,n+1):
      if ni==(N-Npone):
         length=length+1
         end=end+1
      temp=sd(mdata[start:end])
      temparr.append(temp)
      start=end
      end=end+length
   temp=sig(temparr)
   uncertainty.append(temp)
plt.plot(range(1,N+1),uncertainty,'bo-')
plt.show()
print uncertainty
