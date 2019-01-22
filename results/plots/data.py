import os

#Define variables
doscaling=True
doenergies=False
filename1=['he4','o16','snm','ca40']
middle=['lin','ip','quad']
filename2=['.dmc10000','.dmc1000','.dmc1000','.dmc1000_comb']
dos=[3] #which systems to do
dot=[0,1,2] #which types to do
smax=len(dos)
tmax=len(dot)
minval=10

if doscaling:
   print 'Scaling:'
   for ns in range(0,smax): #loop of systems
      s=dos[ns]
      allmean=[]
      for nt in range(0,tmax): #loop of types of correlations
         t=dot[nt]
         filename="../"+filename1[s]+"/"+filename1[s]+middle[t]+filename2[s]+".out"
         os.system("grep 'Time for block' "+filename+" > temp.dat")
         f=open('temp.dat')
         mdata=[]
         for line in f.readlines():
            mdata.append(float(line[38:50]))
         f.close()
         os.system("rm temp.dat")
         allmean.append(sum(mdata)/len(mdata))
      print filename1[s]+middle[1]+'='+str(allmean[1]/allmean[0])
      if tmax==3:
         print filename1[s]+middle[2]+'='+str(allmean[2]/allmean[0])

if doenergies:
   print 'Energies:'
   for ns in range(0,smax): #loop of systems
      s=dos[ns]
      allmean=[]
      for nt in range(0,tmax): #loop of types of correlations
         t=dot[nt]
         filename="../"+filename1[s]+"/"+filename1[s]+middle[t]+filename2[s]+".dat"
         f=open(filename)
         mdata=[]
         edata=[]
         for line in f.readlines():
            mdata.append(float(line[47:64]))
            edata.append(float(line[71:-4]))
         f.close()
         mean=0
         err=0
         count=0
         length=len(mdata)
         for i in range(minval-1,length):
            mean=mean+mdata[i]
            err=err+edata[i]**2
            count=count+1
         mean=mean/(count)
         err=(err/(count))**(0.5)
         print "E"+filename1[s]+middle[t]+" = "+str(round(mean,2))+"+/-"+str(round(err,2))+" MeV"
