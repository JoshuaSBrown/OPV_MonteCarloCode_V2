# Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -pedantic -g

# Targets
ALL : test_linklist

test_linklist : test_linklist.c linklist.o 
	$(CC) $(CFLAGS) -o test_linklist test_linklist.c linklist.o 

linklist.o : linklist.c
	$(CC) $(CFLAGS) -c linklist.c

#mem.o : ../../../MEM/mem.c
#	$(CC) $(CFLAGS) -c ../../../MEM/mem.c

.PHONY: clean
clean :
	$(RM) *.o
