import os

#Define variables
file1= 'o16noip.dmc.1000.dat' #no ind pair
file2= 'dat.6742264' #ind pair
min=1

#plot with gnuplot
os.system('gnuplot -e "file1=\''+file1+'\'; file2=\''+file2+'\'" -persist plot.gp')

#calculations for no ip and print those
mdata=[]
edata=[]
f = open(file1)
for line in f.readlines():
   mdata.append(float(line[47:64]))
   edata.append(float(line[71:-4]))
f.close()

length=len(mdata)
mean=0
err=0
for i in range(min-1,length):
   mean=mean+mdata[i]
   err=err+edata[i]**2
mean=mean/(length)
err=(err/(length))**(0.5)
os.system('echo Enoip = '+str(mean)+' +/- '+str(err)+' MeV')

#calculations for ip and print those
mdata=[]
edata=[]
f = open(file2)
for line in f.readlines():
   mdata.append(float(line[47:64]))
   edata.append(float(line[71:-4]))
f.close()

length=len(mdata)
mean=0
err=0
for i in range(min-1,length):
   mean=mean+mdata[i]
   err=err+edata[i]**2
mean=mean/(length)
err=(err/(length))**(0.5)
os.system('echo Eip = '+str(mean)+' +/- '+str(err)+' MeV')
