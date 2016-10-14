# Variables
CC = gcc
CFLAGS = -lm -Wall -Wextra -std=gnu99 -pedantic -g 

# Targets
ALL : test_clusterfunctions test_clusterfunctions2 test_clusterfunctions3

test_clusterfunctions : test_clusterfunctions.c arbarray.o omll.o electrode.o midpoint.o neighll.o neighnode.o node.o clusterfunctions.o sitenode.o matrix.o linklist.o cluster.o montecarlo.o read.o error.o
	$(CC) $(CFLAGS) test_clusterfunctions.c arbarray.o omll.o  electrode.o midpoint.o neighll.o neighnode.o node.o clusterfunctions.o sitenode.o matrix.o linklist.o cluster.o montecarlo.o read.o error.o -o test_clusterfunctions

test_clusterfunctions2 : test_clusterfunctions2.c arbarray.o omll.o electrode.o midpoint.o neighll.o neighnode.o node.o clusterfunctions.o sitenode.o matrix.o linklist.o cluster.o montecarlo.o read.o error.o
	$(CC) $(CFLAGS) test_clusterfunctions2.c arbarray.o omll.o electrode.o midpoint.o neighll.o neighnode.o node.o clusterfunctions.o sitenode.o matrix.o linklist.o cluster.o montecarlo.o read.o error.o -o test_clusterfunctions2

test_clusterfunctions3 : test_clusterfunctions3.c arbarray.o omll.o electrode.o midpoint.o neighll.o neighnode.o node.o clusterfunctions.o sitenode.o matrix.o linklist.o cluster.o montecarlo.o read.o error.o
	$(CC) $(CFLAGS) test_clusterfunctions3.c midpoint.o arbarray.o omll.o electrode.o neighll.o neighnode.o node.o clusterfunctions.o sitenode.o matrix.o linklist.o cluster.o montecarlo.o read.o error.o -o test_clusterfunctions3

omll.o : ../OMLL/omll.c
	$(CC) $(CFLAGS) -c ../OMLL/omll.c -o omll.o

electrode.o : ../ELECTRODE/electrode.c
	$(CC) $(CFLAGS) -c ../ELECTRODE/electrode.c -o electrode.o

midpoint.o : ../MIDPOINT/midpoint.c
	$(CC) $(CFLAGS) -c ../MIDPOINT/midpoint.c -o midpoint.o

neighnode.o : ../NEIGHNODE/neighnode.c
	$(CC) $(CFLAGS) -c ../NEIGHNODE/neighnode.c -o neighnode.o

node.o : ../NODE/node.c
	$(CC) $(CFLAGS) -c ../NODE/node.c -o node.o

arbarray.o : ../ARBARRAY/arbarray.c
	$(CC) $(CFLAGS) -c ../ARBARRAY/arbarray.c -o arbarray.o

neighll.o : ../NEIGHLL/neighll.c
	$(CC) $(CFLAGS) -c ../NEIGHLL/neighll.c -o neighll.o

montecarlo.o : ../MONTECARLO/montecarlo.c
	$(CC) $(CFLAGS) -c ../MONTECARLO/montecarlo.c -o montecarlo.o

sitenode.o : ../SITENODE/sitenode.c
	$(CC) $(CFLAGS) -c ../SITENODE/sitenode.c -o sitenode.o

matrix.o : ../MATRIX/matrix.c
	$(CC) $(CFLAGS) -c ../MATRIX/matrix.c -o matrix.o

linklist.o : ../LINKLIST/linklist.c
	$(CC) $(CFLAGS) -c ../LINKLIST/linklist.c -o linklist.o

cluster.o : ../CLUSTER/cluster.c
	$(CC) $(CFLAGS) -c ../CLUSTER/cluster.c -o cluster.o

clusterfunctions.o : clusterfunctions.c midpoint.o
	$(CC) $(CFLAGS) midpoint.o -c clusterfunctions.c -o clusterfunctions.o 

read.o : ../PARAMETERS/read.c
	$(CC) $(CFLAGS) -c ../PARAMETERS/read.c -o read.o

error.o : ../ERROR/error.c
	$(CC) $(CFLAGS) -c ../ERROR/error.c -o error.o

.PHONY : clean
clean :
	$(RM) *.o

