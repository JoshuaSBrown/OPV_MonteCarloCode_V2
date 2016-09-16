# Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -pedantic -g 

# Targets
ALL : test_charge

test_charge : test_charge.c charge.o matrix.o matrix_linklist.o sitenode.o 
	$(CC) $(CFLAGS) -o test_charge test_charge.c charge.o matrix.o matrix_linklist.o sitenode.o

charge.o : charge.c
	$(CC) $(CFLAGS) -c charge.c

matrix.o : ../MATRIX/matrix.c
	$(CC) $(CFLAGS) -c ../MATRIX/matrix.c

matrix_linklist.o : ../MATRIX_LINKLIST/matrix_linklist.c
	$(CC) $(CFLAGS) -c ../MATRIX_LINKLIST/matrix_linklist.c

sitenode.o : ../SITENODE/sitenode.c
	$(CC) $(CFLAGS) -c ../SITENODE/sitenode.c

.PHONY: clean
clean:
	$(RM) *.o
