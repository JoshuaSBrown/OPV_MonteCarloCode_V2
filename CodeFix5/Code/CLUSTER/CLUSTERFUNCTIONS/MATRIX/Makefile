# Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -pedantic -g

# Targets
ALL : test_matrix test_matrix2

test_matrix : test_matrix.c matrix.o
	$(CC) $(CFLAGS) -o test_matrix test_matrix.c matrix.o

test_matrix2 : test_matrix2.c matrix.o
	$(CC) $(CFLAGS) -o test_matrix2 test_matrix2.c matrix.o

matrix.o : matrix.c
	$(CC) $(CFLAGS) -c matrix.c

#mem.o : ../../../MEM/mem.c
#	$(CC) $(CFLAGS) -c ../../../MEM/mem.c

.PHONY: clean
clean :
	$(RM) *.o
