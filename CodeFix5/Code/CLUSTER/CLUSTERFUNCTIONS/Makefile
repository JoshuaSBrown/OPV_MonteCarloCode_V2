# Variables
CC = gcc
CFLAGS = -lm -Wall -Wextra -std=gnu99 -pedantic -g 

# Targets
ALL : test_clusterfunctions test_clusterfunctions2 test_clusterfunctions3

test_clusterfunctions : test_clusterfunctions.c clusterfunctions.o sitenode.o matrix.o linklist.o cluster.o montecarlo.o read.o
	$(CC) $(CFLAGS) test_clusterfunctions.c clusterfunctions.o sitenode.o matrix.o linklist.o cluster.o montecarlo.o read.o -o test_clusterfunctions

test_clusterfunctions2 : test_clusterfunctions2.c clusterfunctions.o sitenode.o matrix.o linklist.o cluster.o montecarlo.o read.o
	$(CC) $(CFLAGS) test_clusterfunctions2.c clusterfunctions.o sitenode.o matrix.o linklist.o cluster.o montecarlo.o read.o -o test_clusterfunctions2

test_clusterfunctions3 : test_clusterfunctions3.c clusterfunctions.o sitenode.o matrix.o linklist.o cluster.o montecarlo.o read.o
	$(CC) $(CFLAGS) test_clusterfunctions3.c clusterfunctions.o sitenode.o matrix.o linklist.o cluster.o montecarlo.o read.o -o test_clusterfunctions3

montecarlo.o : MONTECARLO/montecarlo.c
	$(CC) $(CFLAGS) -c MONTECARLO/montecarlo.c -o montecarlo.o

sitenode.o : SITENODE/sitenode.c
	$(CC) $(CFLAGS) -c SITENODE/sitenode.c -o sitenode.o

matrix.o : MATRIX/matrix.c
	$(CC) $(CFLAGS) -c MATRIX/matrix.c -o matrix.o

linklist.o : LINKLIST/linklist.c
	$(CC) $(CFLAGS) -c LINKLIST/linklist.c -o linklist.o

cluster.o : DATASTRUCT/cluster.c
	$(CC) $(CFLAGS) -c DATASTRUCT/cluster.c -o cluster.o

clusterfunctions.o : clusterfunctions.c
	$(CC) $(CFLAGS) -c clusterfunctions.c -o clusterfunctions.o

read.o : ../../PARAMETERS/read.c
	$(CC) $(CFLAGS) -c ../../PARAMETERS/read.c -o read.o

#mem.o : ../../MEM/mem.c
#	$(CC) $(CFLAGS) -c ../../MEM/mem.c -o mem.o

.PHONY : clean
clean :
	$(RM) *.o

