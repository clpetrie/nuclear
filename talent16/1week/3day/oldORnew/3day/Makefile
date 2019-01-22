FC=gfortran
EXEC=qmc
#FLAGS=-O3 #-check all -warn all
FLAGS=-O0 -fbounds-check -ffpe-trap=invalid,zero,overflow -g #-Wall #-check all -warn all
#FLAGS=-O #-check all -warn all

.SUFFIXES: .o .f .f90

.f90.o:
	$(FC) $(FLAGS) -c $<

%.o: %.mod

OBJECTS=\
   ran.o\
   mymod.o\
   qmc.o\

qmc: $(OBJECTS)
	$(FC) -o $(EXEC) $(OBJECTS)
#	$(FC) -o ./$(EXEC) $(OBJECTS)

clean:
	rm -f *\.mod *\.o $(EXEC)
