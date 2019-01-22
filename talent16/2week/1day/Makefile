FC=mpif77
EXEC=qmc
FLAGS=-O3

.SUFFIXES: .o .f .f90

.f90.o:
	$(FC) $(FLAGS) -c $<

%.o: %.mod

OBJECTS=\
   stack.o\
   mympi.o\
   wavefunction.o\
   ran.o\
   estimator.o\
   step.o\
   operators.o\
   qmc.o

qmc: $(OBJECTS)
	$(FC) -o ./$(EXEC) $(OBJECTS)

clean:
	rm -f *\.mod *\.o $(EXEC)
