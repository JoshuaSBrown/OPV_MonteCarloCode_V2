# Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=gnu99 -pedantic -g

# Targets
ALL : test_cluster

test_cluster : test_cluster.c  cluster.o 
	$(CC) $(CFLAGS) test_cluster.c cluster.o -o test_cluster

cluster.o : cluster.c
	$(CC) $(CFLAGS) -c cluster.c -o cluster.o

#mem.o : ../../../MEM/mem.c
#	$(CC) $(CFLAGS) -c ../../../MEM/mem.c -o mem.o

.PHONY : clean
clean :
	$(RM) *.o *exe
