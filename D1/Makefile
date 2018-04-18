# -*- Makefile -*-
SHELL=/bin/sh
CC=gfortran
CFLAGS=-ggdb -Wall -std=f2003


default:utilities.o histogram.x montecarlo.x

utilities.o: utilities.f03
	$(CC) -o $@ -c $< $(CFLAGS)

montecarlo.o: montecarlo.f03
	$(CC) -o $@ -c $< $(CFLAGS)

montecarlo.x: montecarlo.o utilities.o
	$(CC) -o $@ $^ $(CFLAGS)

histogram.o: histogram.f03
	$(CC) -o $@ -c $< $(CFLAGS)

histogram.x: histogram.o utilities.o
	$(CC) -o $@ $^ $(CFLAGS)


clean:
	rm -f comp_* *.mod ./src/*.mod ./src/*.o *.o *.x *.out /src/*.out ./src/*~ *~
plot:
	gnuplot plot.gp

check: ./run.sh
	./run.sh
	cmp ./reference/sum_int comp_sum_int || exit 1
	rm comp_sum_int
	cmp ./reference/sum_real comp_sum_real || exit 1
	rm comp_sum_real
	cmp ./reference/sort_int comp_sort_int || exit 1
	rm comp_sort_int
	cmp ./reference/sort_real comp_sort_real || exit 1
	rm comp_sort_real
	cmp comp_sum_derived ./reference/sum_derived || exit 1
	rm comp_sum_derived
