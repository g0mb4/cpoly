all:
	gcc -o polytest test.c poly.c -lm -lgsl -lgslcblas
