sho:mymod.o sho.o
	gfortran -o sho sho.o mymod.o

sho.o: sho.f90
	gfortran -O -c sho.f90

mymod.o: mymod.f90
	gfortran -O -c mymod.f90
