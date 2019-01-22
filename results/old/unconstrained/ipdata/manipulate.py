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
   filename="he4ip.unc10000_"+str(i)+".out"
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
uin=raw_input('Please enter a max value between '+str(imin)+' and '+str(tot)+': ')
reject=True
while reject:
   if str(uin).isdigit() and int(uin)>=imin and int(uin)<=tot:
      imax=int(uin)
      reject=False
   else:
      print 'Please enter an integer!'
      uin=raw_input('Please enter a max value between '+str(imin)+' and '+str(tot)+': ')

while True:
   plt.clf()
   avedata=[]
   errdata=[]
   length=len(data[:][0])
   if tdata==[]:
      tdata=range(1,length+1)
   for i in range(0,length):
      total=[]
      for j in range(imin-1,imax):
         total.append(data[j][i])
      avedata.append(average(total))
      errdata.append(sig(total))
#      avedata.append(np.mean(total))
#      errdata.append(np.std(total))
#   print avedata
#   print errdata
   plt.errorbar(tdata,avedata,yerr=errdata,fmt='o')
   plt.xlabel(r'$\tau$')
   plt.ylabel(r'E')
   plt.ylim([-30,-24])
   plt.show(block=False)

   print 'Imax = '+str(imax)
   uin=raw_input('Input: ')
   redo=True
   while redo:
      if uin=='help':
         print "Please enter 'p' (plus), 'm' (minus), 'exit', or a number for imax between "+str(imin)+" and "+str(tot)+": "
         uin=raw_input('Input: ')
      elif uin=='p':
         if imax<tot:
            imax=imax+1
            redo=False
         else:
            print "Already at max"
            uin=raw_input('Input: ')
      elif uin=='m':
         if imax>imin:
            imax=imax-1
            redo=False
         else:
            print "Already at min"
            uin=raw_input('Input: ')
      elif uin=='exit':
         sys.exit()
      elif str(uin).isdigit():
         imax=int(uin)
         redo=False
      else:
         print "Please enter 'p' (plus), 'm' (minus), 'exit', or a number for imax between "+str(imin)+" and "+str(tot)+": "
         uin=raw_input('Input: ')
