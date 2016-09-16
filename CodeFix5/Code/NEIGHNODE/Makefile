# Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=gnu99 -pedantic -g

# Targets
ALL : test_neighnode

test_neighnode : test_neighnode.c  neighnode.o 
	$(CC) $(CFLAGS) test_neighnode.c neighnode.o -o test_neighnode

neighnode.o : neighnode.c
	$(CC) $(CFLAGS) -c neighnode.c -o neighnode.o

.PHONY : clean
clean :
	$(RM) *.o *exe
