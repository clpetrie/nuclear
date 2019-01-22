import os

filename1="../ca40/ca40lin.dmc1000_comb.out"
filename2="../ca40/ca40ip.dmc1000_comb.out"
#filename2="../ca40/ca40quad.dmc1000_comb.out"
os.system("grep 'Time for block' "+filename1+" > temp1.dat")
os.system("grep 'Time for block' "+filename2+" > temp2.dat")
f=open('temp1.dat')
mdata1=[]
mdata2=[]
for line in f.readlines():
   mdata1.append(float(line[42:50]))
f.close()
f=open('temp2.dat')
mdata=[]
for line in f.readlines():
   mdata2.append(float(line[38:50]))
f.close()
#os.system("rm temp1.dat")
#os.system("rm temp2.dat")
print mdata1
print ''
print ''
print ''
print mdata2
mean1=sum(mdata1)/len(mdata1)
mean2=sum(mdata2)/len(mdata2)
print mean2/mean1
