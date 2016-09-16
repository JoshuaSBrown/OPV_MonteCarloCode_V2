# Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=gnu99 -pedantic -g

# Targets
ALL : test_omll

test_omll : test_omll.c  omll.o midpoint.o
	$(CC) $(CFLAGS) test_omll.c omll.o midpoint.o -o test_omll

omll.o : omll.c
	$(CC) $(CFLAGS) -c omll.c -o omll.o

midpoint.o : ../MIDPOINT/midpoint.c
	$(CC) $(CFLAGS) -c ../MIDPOINT/midpoint.c -o midpoint.o


.PHONY : clean
clean :
	$(RM) *.o *exe
