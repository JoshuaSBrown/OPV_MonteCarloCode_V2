# Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -pedantic -g

# Targets
ALL : test_linklist2

test_linklist2 : test_linklist2.c linklist2.o 
	$(CC) $(CFLAGS) -o test_linklist2 test_linklist2.c linklist2.o 

linklist2.o : linklist2.c
	$(CC) $(CFLAGS) -c linklist2.c

.PHONY: clean
clean :
	$(RM) *.o
