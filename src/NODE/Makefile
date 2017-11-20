# Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=gnu99 -pedantic -g

# Targets
ALL : test_node

test_node : test_node.c  node.o 
	$(CC) $(CFLAGS) test_node.c node.o -o test_node

node.o : node.c
	$(CC) $(CFLAGS) -c node.c -o node.o

.PHONY : clean
clean :
	$(RM) *.o *exe
