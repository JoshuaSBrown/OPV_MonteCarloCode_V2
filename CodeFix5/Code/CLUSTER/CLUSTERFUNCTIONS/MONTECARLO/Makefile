# Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=gnu99 -pedantic -g
LDLIBS = -lm

# Targets
ALL : test_montecarlo

test_montecarlo : test_montecarlo.c montecarlo.o sitenode.o matrix.o 
	$(CC) $(CFLAGS) $(LDLIBS) -o test_montecarlo test_montecarlo.c montecarlo.o sitenode.o matrix.o 

montecarlo.o : montecarlo.c 
	$(CC) $(CFLAGS) $(LDLIBS) -c montecarlo.c -o montecarlo.o

sitenode.o : ../SITENODE/sitenode.c 
	$(CC) $(CFLAGS) -c  ../SITENODE/sitenode.c -o sitenode.o 

matrix.o : ../MATRIX/matrix.c 
	$(CC) $(CFLAGS)  -c ../MATRIX/matrix.c -o matrix.o

#mem.o : ../../../MEM/mem.c
#	$(CC) $(CFLAGS) -c ../../../MEM/mem.c -o mem.o

.PHONY: clean
clean:
	$(RM) *.o 

