MPI=true

ifeq ($(MPI),true)
   FC=mpif90
else
   FC=gfortran
endif

EXEC=day4
FLAGS=-O3 #-check all -warn all
#FLAGS=-O0 -fbounds-check -ffpe-trap=invalid,zero,overflow -g -fcheck=all #-Wall #-check all -warn all
#FLAGS=-O #-check all -warn all

.SUFFIXES: .o .f .f90

.f90.o:
#	$(FC) $(FLAGS) -c $<
	$(FC) $(FLAGS) -c -I/usr/local/include $<

%.o: %.mod

OBJECTS=\
   ran.o\
   step.o\
   day4.o\

day4: $(OBJECTS)
#	$(FC) -o $(EXEC) $(OBJECTS)
	$(FC) -o ./$(EXEC) $(OBJECTS)

clean:
	rm -f *\.mod *\.o $(EXEC)
