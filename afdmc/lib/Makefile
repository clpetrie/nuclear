#$Id: Makefile,v 1.2 2003/06/26 13:21:06 nuclear Exp $
#Intel fortran compiler -xW requires pentium 4 instructions
#Use -tpp6 for pentium 3
FC=ifort
FFLAGS=-O3

#FC=gfortran
#FFLAGS=-O3 #-march=pentium4 -mtune=pentium4

NOOPT=-O0
AR=ar
RANLIB=ranlib

all:
	(cd blas;\
        make FC="$(FC)" FFLAGS="$(FFLAGS)" AR="$(AR)" RANLIB="$RANLIB" )
	(cd lapack;\
        make FC="$(FC)" FFLAGS="$(FFLAGS)" AR="$(AR)" RANLIB="$RANLIB"\
            NOOPT="$(NOOPT)" )
	cp blas/libblasmin.a .
	cp lapack/liblapackmin.a .

clean:
	(cd blas; make clean)
	(cd lapack; make clean)

veryclean:
	rm -rf liblapackmin.a libblasmin.a
	(cd blas; make clean)
	(cd lapack; make clean)
	
