import numpy as np
import matplotlib.pyplot as plt
num_iter=100
pi=3.14159
radius=1.0
Abox=4*radius**2
points_list=[]
mean_list=[]
std_list=[]
minnum=10
iternum=minnum
maxnum=1000+iternum

for num in range(minnum,maxnum,iternum):
   pi_list=[]
   for i in range(0,num_iter):
      inside=0
      for j in range(0,num):
         x,y=2*radius*np.random.random(2)-radius
         r=np.sqrt(x**2+y**2)
         if r <= radius:
            inside=inside+1

      pi_list.append(float(inside)/float(num)*Abox)
   points_list.append(num)
   mean_list.append(np.mean(pi_list))
   std_list.append(np.std(pi_list))
   print "Completed num = "+str(num)+"/"+str(maxnum-iternum)

f = plt.figure()
plt.xlabel("Sampled Points",fontsize=20)
plt.ylabel(r"Estimate of $\pi$",fontsize=20)
plt.xlim([min(points_list)-10, max(points_list)+10])
plt.axhline(linewidth=3, color='r',y=pi)
plt.errorbar(points_list,mean_list,yerr=std_list,fmt='o')
plt.show()

#f.savefig("pi_estimate.pdf",bbox_inches='tight')
