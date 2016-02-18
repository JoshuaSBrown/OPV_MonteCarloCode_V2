# Variables
CC = gcc
CFLAGS = -lm -Wall -Wextra -std=gnu99 -pedantic -g 

# Targets
ALL : test_clustersitenode

test_clustersitenode : test_clustersitenode.c clustersitenode.o cluster.o sitenode.o clusterfunctions.o read.o matrix.o linklist.o montecarlo.o
	$(CC) $(CFLAGS) test_clustersitenode.c clustersitenode.o cluster.o sitenode.o clusterfunctions.o read.o matrix.o linklist.o montecarlo.o -o test_clustersitenode

clustersitenode.o : clustersitenode.c
	$(CC) $(CFLAGS) -c clustersitenode.c -o clustersitenode.o

cluster.o : ../CLUSTERFUNCTIONS/DATASTRUCT/cluster.c
	$(CC) $(CFLAGS) -c ../CLUSTERFUNCTIONS/DATASTRUCT/cluster.c -o cluster.o

sitenode.o : ../CLUSTERFUNCTIONS/SITENODE/sitenode.c
	$(CC) $(CFLAGS) -c ../CLUSTERFUNCTIONS/SITENODE/sitenode.c -o sitenode.o

#mem.o : ../../MEM/mem.c
#	$(CC) $(CFLAGS) -c ../../MEM/mem.c -o mem.o

clusterfunctions.o : ../CLUSTERFUNCTIONS/clusterfunctions.c
	$(CC) $(CFLAGS) -c ../CLUSTERFUNCTIONS/clusterfunctions.c -o clusterfunctions.o

matrix.o : ../CLUSTERFUNCTIONS/MATRIX/matrix.c
	$(CC) $(CFLAGS) -c ../CLUSTERFUNCTIONS/MATRIX/matrix.c -o matrix.o

linklist.o : ../CLUSTERFUNCTIONS/LINKLIST/linklist.c
	$(CC) $(CFLAGS) -c ../CLUSTERFUNCTIONS/LINKLIST/linklist.c -o linklist.o

montecarlo.o : ../CLUSTERFUNCTIONS/MONTECARLO/montecarlo.c
	$(CC) $(CFLAGS) -c ../CLUSTERFUNCTIONS/MONTECARLO/montecarlo.c -o montecarlo.o

read.o : ../../PARAMETERS/read.c
	$(CC) $(CFLAGS) -c ../../PARAMETERS/read.c -o read.o

.PHONY : clean
clean :
	$(RM) *.o
