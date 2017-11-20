# Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=gnu99 -pedantic -g

# Targets
ALL : test_midpoint

test_midpoint : test_midpoint.c  midpoint.o 
	$(CC) $(CFLAGS) test_midpoint.c midpoint.o -o test_midpoint

midpoint.o : midpoint.c
	$(CC) $(CFLAGS) -c midpoint.c -o midpoint.o

#mem.o : ../../../MEM/mem.c
#	$(CC) $(CFLAGS) -c ../../../MEM/mem.c -o mem.o

.PHONY : clean
clean :
	$(RM) *.o *exe
