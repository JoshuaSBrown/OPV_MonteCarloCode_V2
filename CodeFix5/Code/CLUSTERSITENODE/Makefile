# Variables
CC = gcc
CFLAGS = -lm -Wall -Wextra -std=gnu99 -pedantic -g 

# Targets
ALL : test_clustersitenode test_clustersitenode2

test_clustersitenode : test_clustersitenode.c clustersitenode.o omll.o neighll.o neighnode.o midpoint.o electrode.o read.o node.o arbarray.o cluster.o sitenode.o clusterfunctions.o matrix.o linklist.o montecarlo.o
	$(CC) $(CFLAGS) test_clustersitenode.c clustersitenode.o omll.o neighll.o neighnode.o midpoint.o electrode.o read.o node.o arbarray.o cluster.o sitenode.o clusterfunctions.o matrix.o linklist.o montecarlo.o -o test_clustersitenode

test_clustersitenode2 : test_clustersitenode2.c clustersitenode.o omll.o neighll.o neighnode.o midpoint.o electrode.o read.o node.o arbarray.o cluster.o sitenode.o clusterfunctions.o matrix.o linklist.o montecarlo.o
	$(CC) $(CFLAGS) test_clustersitenode2.c clustersitenode.o omll.o neighll.o neighnode.o midpoint.o electrode.o read.o node.o arbarray.o cluster.o sitenode.o clusterfunctions.o matrix.o linklist.o montecarlo.o -o test_clustersitenode2

clustersitenode.o : clustersitenode.c
	$(CC) $(CFLAGS) -c clustersitenode.c -o clustersitenode.o

cluster.o : ../CLUSTER/cluster.c
	$(CC) $(CFLAGS) -c ../CLUSTER/cluster.c -o cluster.o

node.o : ../NODE/node.c
	$(CC) $(CFLAGS) -c ../NODE/node.c -o node.o

neighnode.o : ../NEIGHNODE/neighnode.c
	$(CC) $(CFLAGS) -c ../NEIGHNODE/neighnode.c -o neighnode.o

neighll.o : ../NEIGHLL/neighll.c
	$(CC) $(CFLAGS) -c ../NEIGHLL/neighll.c -o neighll.o

midpoint.o : ../MIDPOINT/midpoint.c
	$(CC) $(CFLAGS) -c ../MIDPOINT/midpoint.c -o midpoint.o

electrode.o : ../ELECTRODE/electrode.c
	$(CC) $(CFLAGS) -c ../ELECTRODE/electrode.c -o electrode.o

sitenode.o : ../SITENODE/sitenode.c
	$(CC) $(CFLAGS) -c ../SITENODE/sitenode.c -o sitenode.o

clusterfunctions.o : ../CLUSTERFUNCTIONS/clusterfunctions.c
	$(CC) $(CFLAGS) -c ../CLUSTERFUNCTIONS/clusterfunctions.c -o clusterfunctions.o

matrix.o : ../MATRIX/matrix.c
	$(CC) $(CFLAGS) -c ../MATRIX/matrix.c -o matrix.o

linklist.o : ../LINKLIST/linklist.c
	$(CC) $(CFLAGS) -c ../LINKLIST/linklist.c -o linklist.o

montecarlo.o : ../MONTECARLO/montecarlo.c
	$(CC) $(CFLAGS) -c ../MONTECARLO/montecarlo.c -o montecarlo.o

arbarray.o : ../ARBARRAY/arbarray.c
	$(CC) $(CFLAGS) -c ../ARBARRAY/arbarray.c -o arbarray.o

omll.o : ../OMLL/omll.c
	$(CC) $(CFLAGS) -c ../OMLL/omll.c -o omll.o

read.o : ../PARAMETERS/read.c
	$(CC) $(CFLAGS) -c ../PARAMETERS/read.c -o read.o

.PHONY : clean
clean :
	$(RM) *.o
