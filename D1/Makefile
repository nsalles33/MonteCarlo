# -*- Makefile -*-
# SHELL=/bin/sh
# CC=gfortran
# CFLAGS=-ggdb -Wall -std=f2003 -fbounds-check -fcheck=all -I/usr/include
# LD_FLAGS = -L/usr/lib

SHELL=/bin/sh
CC=gfortran
CFLAGS=-ggdb -Wall -std=f2003 -I/usr/include
LD_FLAGS = -L/usr/lib 

default:utilities.o histogram.x

utilities.o: utilities.f03
	$(CC) -o $@ -c $< $(CFLAGS)

histogram.o: histogram.f03
	$(CC) -o $@ -c $< $(CFLAGS)

histogram.x: histogram.o utilities.o
	$(CC) -o $@ $^ $(CFLAGS)


clean:
	rm -f comp_* *.mod ./src/*.mod ./src/*.o *.o *.x *.out /src/*.out ./src/*~ *~
