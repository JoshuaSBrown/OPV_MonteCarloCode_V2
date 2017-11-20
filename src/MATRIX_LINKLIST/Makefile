# Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -pedantic -g

# Targets
ALL : test_matrix_linklist

test_matrix_linklist : test_matrix_linklist.c matrix_linklist.o matrix.o
	$(CC) $(CFLAGS) -o test_matrix_linklist test_matrix_linklist.c matrix_linklist.o matrix.o

matrix_linklist.o : matrix_linklist.c
	$(CC) $(CFLAGS) -c matrix_linklist.c

matrix.o : ..//MATRIX/matrix.c
	$(CC) $(CFLAGS) -c ../MATRIX/matrix.c

.PHONY: clean
clean :
	$(RM) *.o
