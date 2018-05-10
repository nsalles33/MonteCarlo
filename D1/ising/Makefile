#Dirs
EXEDIR=./
SRCDIR=src
HEADDIR=include

MAIN=montecarlo
HEAD=utilities
#Compiler
#OPENBLAS_LIB=/usr/lib/x86_64-linux-gnu/openblas
LFLAGS= -lm -lopenblas -L ${OPENBLAS_LIB}
OPENMP_INC=/u/shared/programs/x86_64/gcc/4.8.2/lib/gcc/x86_64-redhat-linux-gnu/4.8.2/include
CC=gcc
CFLAGS=-g -std=c11 -fopenmp -D_USE_OMP -Wall -Wextra -I include -I ${OPENBLAS_INC} #-I ${OPENMP_INC}

default: montecarlo.x

run: montecarlo.x
	./montecarlo.x

montecarlo.x: ${MAIN}.o ${SRCDIR}/${HEAD}.o
	$(CC) $(CFLAGS) $^ -o ${EXEDIR}/montecarlo.x $(LFLAGS)

${SRCDIR}/${HEAD}.o: ${HEADDIR}/${HEAD}.h


%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@ $(LFLAGS)

clean:
	rm -rf ${SRCDIR}/*.o  *.o *~ ./*.x
