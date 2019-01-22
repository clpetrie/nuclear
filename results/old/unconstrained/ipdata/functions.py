def square(list):
   slist=[]
   for i in list:
      slist.append(i**2)
   return slist

def average(indata):
   ave=0
   length=len(indata)
   for i in range(0,length):
      ave=ave+indata[i]
   ave=ave/length
   return ave

def sig(indata):
   sigma=0
   length=len(indata)
   mean=average(indata)
   mean2=average(square(indata))
   sigma=(1.0/length*abs(mean2-mean**2))**(0.5)
   return sigma
