day1:mymod.o day1.o
	gfortran -o day1 day1.o mymod.o

day1.o: day1.f90
	gfortran -O -c day1.f90

mymod.o: mymod.f90
	gfortran -O -c mymod.f90
