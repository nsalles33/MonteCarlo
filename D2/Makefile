# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -std=c11 -DMPI -I ../include
#LDLIBS=-lm

# list of source files
SRC=lcg_rand_gen.c
VPATH=.
OBJ=lcg_rand_gen.o #input.o main.o output.o utilities.o verlet1.o verlet2.o

default: lcg_rand_gen.x

clean:
	rm -f *.mod *.o ../ljmd-mpi.x

# linker rule
lcg_rand_gen.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
