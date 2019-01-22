import os
import matplotlib.pyplot as plt

#Define variables
filename=["o16scaling_w10000_proc40.out", "o16scaling_w10000_proc160.out", "o16scaling_w10000_proc1280.out", "o16scaling_w10000_proc2560.out"]
#filename=["out.6809063", "out.6809070","out.6809078", "out.6809081"]
#filename=["2out.6765308", "3out.6767487", "4out.6767489", "5out.6767494"]
nproc=[40,160,1280,2560]
#nproc=[16,128,1024,2048]

avei=[]
for s in range(0,len(filename)):
   data=[]
   os.system("grep 'Time for propagation' "+str(filename[s])+" > temp.out")
   f = open("temp.out")
   for line in f.readlines():
      data.append(float(line[40:50]))
   f.close()
   os.system("rm temp.out")
   ave=sum(data)/float(len(data))
   avei.append(1/ave)

   print str(nproc[s])+", "+str(1/ave)
#   print str(nproc[s])+", "+str(ave)+", "+str(nproc[s]*ave)

plt.xlabel('# of cores')
plt.ylabel(r'1/t (sec$^{-1}$)')
plt.plot(nproc,avei,marker='o')
plt.show()
