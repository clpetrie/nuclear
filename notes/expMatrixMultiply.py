from numpy import *
from numpy.polynomial.hermite import hermroots, hermval
from math import factorial
from cmath import sqrt
from matplotlib.pyplot import plot, show

N=3
tot=zeros(N)
print tot
print tot[1]
print ''
print ''
print 'blah'
for Ni in range(0,N):
   a=0.01j
   O=5.34
   pi=3.14159
   i=1j
   coef=zeros(N+1)
   coefm1=zeros(N+1)
   for n in range(0,N+1):
      if n==N-1:
          coefm1[n]=1
      if n==N:
          coef[n]=1

   roots = hermroots(coef)
   wi=2**(N-1)*factorial(N)*sqrt(pi)/(N**2*(hermval(roots,coefm1))**2)

   print 'n=',N
   print 'roots=',roots
   print 'wi=',wi

   print ''
   print exp(-a*O**2)

   #solve using Hubbard-Stratanovich and Hermite-Gauss Quadrature
   for n in range(0,N):
      tot[Ni]=tot[Ni]+wi[n]*exp(2*sqrt(-a)*roots[n]*O)
   tot[Ni] = (1.0/pi)**0.5*tot[Ni]

print tot
#plot(tot)
#show()
