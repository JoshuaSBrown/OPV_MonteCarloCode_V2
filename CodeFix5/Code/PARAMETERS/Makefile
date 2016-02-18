# Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -pedantic -g 

# Targets
ALL : test_read

test_read : test_read.c read.o
	$(CC) $(CFLAGS) -o test_read test_read.c read.o 

read.o : read.c
	$(CC) $(CFLAGS) -c read.c

#mem.o : ../MEM/mem.c
#	$(CC) $(CFLAGS) -c ../MEM/mem.c

.PHONY: clean
clean:
	$(RM) *.o
