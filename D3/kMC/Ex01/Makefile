SHELL=/bin/sh

CC=gcc
F90=gfortran

FFLAGS= -Wall -cpp -shared -DSHARED_LIB -fPIC
CFLAGS= -fPIC -shared

SRC=../src

all: fort c run

fort:
	$(F90) $(FFLAGS) -o lib_Ex01.so $(SRC)/lev00_new_types.f90 $(SRC)/lev00_module.f90 $(SRC)/lev00_tools.f90 Ex01.f90

c:
	$(F90) -c -DSHARED_LIB -fPIC $(SRC)/lev00_new_types.f90
	$(CC) $(CFLAGS) -o lib_Ex01c.so lev00_new_types.o Ex01.c

run:
	cp $(SRC)/KMC.x .
	./KMC.x input_KMC.dat

clean:
	rm *.o *.mod *.so *.x
