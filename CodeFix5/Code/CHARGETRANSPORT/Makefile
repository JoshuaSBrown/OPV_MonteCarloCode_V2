# Variables
CC = gcc
CFLAGS = -lm -Wall -Wextra -std=gnu99 -pedantic -g 

# Targets
ALL : test_chargetransport test_chargetransport2

test_chargetransport : test_chargetransport.c chargetransport.o io.o omll.o matrix_linklist.o midpoint.o node.o neighll.o neighnode.o arbarray.o electrode.o matrix.o charge.o linklist.o clusterfunctions.o clustersitenode.o cluster.o sitenode.o montecarlo.o functions.o read.o
	$(CC) $(CFLAGS) test_chargetransport.c chargetransport.o io.o omll.o matrix_linklist.o midpoint.o node.o neighll.o neighnode.o arbarray.o electrode.o charge.o matrix.o linklist.o clusterfunctions.o clustersitenode.o cluster.o sitenode.o functions.o montecarlo.o read.o -o test_chargetransport

test_chargetransport2 : test_chargetransport2.c chargetransport.o io.o omll.o matrix_linklist.o midpoint.o node.o neighll.o neighnode.o arbarray.o electrode.o matrix.o charge.o linklist.o clusterfunctions.o clustersitenode.o cluster.o sitenode.o montecarlo.o functions.o read.o
	$(CC) $(CFLAGS) test_chargetransport2.c chargetransport.o io.o omll.o matrix_linklist.o midpoint.o node.o neighll.o neighnode.o arbarray.o electrode.o charge.o matrix.o linklist.o clusterfunctions.o clustersitenode.o cluster.o sitenode.o functions.o montecarlo.o read.o -o test_chargetransport2

chargetransport.o  : chargetransport.c
	$(CC) $(CFLAGS) -c chargetransport.c -o chargetransport.o

charge.o 					 : ../CHARGE/charge.c
	$(CC) $(CFLAGS) -c ../CHARGE/charge.c -o charge.o

functions.o 			 : ../FUNCTIONS/functions.c
	$(CC) $(CFLAGS) -c ../FUNCTIONS/functions.c -o functions.o

clusterfunctions.o : ../CLUSTERFUNCTIONS/clusterfunctions.c
	$(CC) $(CFLAGS) -c ../CLUSTERFUNCTIONS/clusterfunctions.c -o clusterfunctions.o

linklist.o 				 : ../LINKLIST/linklist.c
	$(CC) $(CFLAGS) -c ../LINKLIST/linklist.c -o linklist.o

montecarlo.o 			 : ../MONTECARLO/montecarlo.c
	$(CC) $(CFLAGS) -c ../MONTECARLO/montecarlo.c -o montecarlo.o

sitenode.o 				 : ../SITENODE/sitenode.c
	$(CC) $(CFLAGS) -c ../SITENODE/sitenode.c -o sitenode.o

cluster.o 				 : ../CLUSTER/cluster.c
	$(CC) $(CFLAGS) -c ../CLUSTER/cluster.c -o cluster.o

node.o 						 : ../NODE/node.c
	$(CC) $(CFLAGS) -c ../NODE/node.c -o node.o

neighll.o 				 : ../NEIGHLL/neighll.c
	$(CC) $(CFLAGS) -c ../NEIGHLL/neighll.c -o neighll.o

neighnode.o 			 : ../NEIGHNODE/neighnode.c
	$(CC) $(CFLAGS) -c ../NEIGHNODE/neighnode.c -o neighnode.o

omll.o 						 : ../OMLL/omll.c
	$(CC) $(CFLAGS) -c ../OMLL/omll.c -o omll.o

electrode.o 	     : ../ELECTRODE/electrode.c
	$(CC) $(CFLAGS) -c ../ELECTRODE/electrode.c -o electrode.o

arbarray.o 				 : ../ARBARRAY/arbarray.c
	$(CC) $(CFLAGS) -c ../ARBARRAY/arbarray.c -o arbarray.o

matrix.o 					 : ../MATRIX/matrix.c
	$(CC) $(CFLAGS) -c ../MATRIX/matrix.c -o matrix.o

clustersitenode.o  : ../CLUSTERSITENODE/clustersitenode.c
	$(CC) $(CFLAGS) -c ../CLUSTERSITENODE/clustersitenode.c -o clustersitenode.o

midpoint.o 				 : ../MIDPOINT/midpoint.c
	$(CC) $(CFLAGS) -c ../MIDPOINT/midpoint.c -o midpoint.o

read.o 						 : ../PARAMETERS/read.c
	$(CC) $(CFLAGS) -c ../PARAMETERS/read.c -o read.o

matrix_linklist.o  : ../MATRIX_LINKLIST/matrix_linklist.c
	$(CC) $(CFLAGS) -c ../MATRIX_LINKLIST/matrix_linklist.c -o matrix_linklist.o

io.o 							 : ../IO/io.c
	$(CC) $(CFLAGS) -c ../IO/io.c -o io.o

.PHONY: clean
clean:
	$(RM) *.o *.exe
