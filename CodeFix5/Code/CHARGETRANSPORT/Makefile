# Variables
CC = gcc
CFLAGS = -lm -Wall -Wextra -std=gnu99 -pedantic -g 

# Targets
ALL : test_chargetransport test_chargetransport2

test_chargetransport : test_chargetransport.c chargetransport.o matrix.o charge.o linklist.o clusterfunctions.o clustersitenode.o cluster.o sitenode.o montecarlo.o functions.o read.o
	$(CC) $(CFLAGS) test_chargetransport.c chargetransport.o charge.o matrix.o linklist.o clusterfunctions.o clustersitenode.o cluster.o sitenode.o functions.o montecarlo.o read.o -o test_chargetransport

test_chargetransport2 : test_chargetransport2.c chargetransport.o matrix.o charge.o linklist.o clusterfunctions.o clustersitenode.o cluster.o sitenode.o montecarlo.o functions.o read.o
	$(CC) $(CFLAGS) test_chargetransport2.c chargetransport.o charge.o matrix.o linklist.o clusterfunctions.o clustersitenode.o cluster.o sitenode.o functions.o montecarlo.o read.o -o test_chargetransport2

chargetransport.o : chargetransport.c
	$(CC) $(CFLAGS) -c chargetransport.c -o chargetransport.o

charge.o : ../CHARGE/charge.c
	$(CC) $(CFLAGS) -c ../CHARGE/charge.c -o charge.o

#mem.o : ../MEM/mem.c
#	$(CC) $(CFLAGS) -c ../MEM/mem.c -o mem.o

functions.o : ../FUNCTIONS/functions.c
	$(CC) $(CFLAGS) -c ../FUNCTIONS/functions.c -o functions.o

clusterfunctions.o : ../CLUSTER/CLUSTERFUNCTIONS/clusterfunctions.c
	$(CC) $(CFLAGS) -c ../CLUSTER/CLUSTERFUNCTIONS/clusterfunctions.c -o clusterfunctions.o

linklist.o : ../CLUSTER/CLUSTERFUNCTIONS/LINKLIST/linklist.c
	$(CC) $(CFLAGS) -c ../CLUSTER/CLUSTERFUNCTIONS/LINKLIST/linklist.c -o linklist.o

montecarlo.o : ../CLUSTER/CLUSTERFUNCTIONS/MONTECARLO/montecarlo.c
	$(CC) $(CFLAGS) -c ../CLUSTER/CLUSTERFUNCTIONS/MONTECARLO/montecarlo.c -o montecarlo.o

sitenode.o : ../CLUSTER/CLUSTERFUNCTIONS/SITENODE/sitenode.c
	$(CC) $(CFLAGS) -c ../CLUSTER/CLUSTERFUNCTIONS/SITENODE/sitenode.c -o sitenode.o

cluster.o : ../CLUSTER/CLUSTERFUNCTIONS/DATASTRUCT/cluster.c
	$(CC) $(CFLAGS) -c ../CLUSTER/CLUSTERFUNCTIONS/DATASTRUCT/cluster.c -o cluster.o

matrix.o : ../CLUSTER/CLUSTERFUNCTIONS/MATRIX/matrix.c
	$(CC) $(CFLAGS) -c ../CLUSTER/CLUSTERFUNCTIONS/MATRIX/matrix.c -o matrix.o

clustersitenode.o : ../CLUSTER/CLUSTERSITENODE/clustersitenode.c
	$(CC) $(CFLAGS) -c ../CLUSTER/CLUSTERSITENODE/clustersitenode.c -o clustersitenode.o

read.o : ../PARAMETERS/read.c
	$(CC) $(CFLAGS) -c ../PARAMETERS/read.c -o read.o

.PHONY: clean
clean:
	$(RM) *.o
