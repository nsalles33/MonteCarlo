SHELL=/bin/sh

CC=
F90=

FFLAGS= -Wall -cpp -shared -DSHARED_LIB -fPIC 
CFLAGS= -fPIC -shared

SRC=../src

all: fort c run

fort:
	$(F90) $(FFLAGS) -o lib_Ex02.so $(SRC)/lev00_new_types.f90 $(SRC)/lev00_module.f90 $(SRC)/lev00_tools.f90 Ex02.f90

c:
	$(F90) -c -DSHARED_LIB -fPIC $(SRC)/lev00_new_types.f90 
	$(CC) $(CFLAGS) -o lib_Ex02c.so lev00_new_types.o Ex02.c

run:
	cp $(SRC)/KMC.x .
	./KMC.x input_KMC.dat

clean:
	rm *.o *.mod *.so


