# -*- Makefile -*-
SHELL=/bin/sh
CC=gfortran
CFLAGS=-ggdb -Wall -std=f2003


default:list_types.o list_tools.o 01_readint.x 02_readfloat.x 03_is_sorted.x 04_is_sorted_real.x 07_derived_type.x

list_types.o: ./list_types.f03
	$(CC) -c $(CFLAGS) -o $@ list_types.f03

list_tools.o: list_tools.f03
	$(CC) -c $(CFLAGS) -o $@ list_tools.f03

01_readint.o: 01_readint.f03
	$(CC) -o $@ -c $< $(CFLAGS)

01_readint.x: 01_readint.o
	$(CC) -o $@ $< $(CFLAGS)

clean:
	rm -f comp_* *.mod ./src/*.mod ./src/*.o *.o *.x *.out /src/*.out ./src/*~ *~

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
