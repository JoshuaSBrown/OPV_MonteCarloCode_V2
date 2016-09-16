# Variables
CC = gcc
CFLAGS = -lm  -m64 -o0 -Wall -Wextra -std=gnu99 -pedantic -ggdb 

# Targets
ALL : run

run : main.c chargetransport.o midpoint.o electrode.o matrix.o read.o charge.o matrix_linklist.o neighll.o neighnode.o node.o omll.o arbarray.o linklist.o clusterfunctions.o clustersitenode.o cluster.o sitenode.o montecarlo.o functions.o io.o
	$(CC) $(CFLAGS) main.c chargetransport.o midpoint.o electrode.o read.o charge.o matrix_linklist.o neighll.o neighnode.o node.o omll.o arbarray.o matrix.o linklist.o clusterfunctions.o clustersitenode.o cluster.o sitenode.o functions.o io.o montecarlo.o -o run

chargetransport.o  : ../CHARGETRANSPORT/chargetransport.c
	$(CC) $(CFLAGS) -c ../CHARGETRANSPORT/chargetransport.c -o chargetransport.o

charge.o 					 : ../CHARGE/charge.c
	$(CC) $(CFLAGS) -c ../CHARGE/charge.c -o charge.o

functions.o 			 : ../FUNCTIONS/functions.c
	$(CC) $(CFLAGS) -c ../FUNCTIONS/functions.c -o functions.o

read.o 						 : ../PARAMETERS/read.c
	$(CC) $(CFLAGS) -c ../PARAMETERS/read.c -o read.o

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

midpoint.o  			 : ../MIDPOINT/midpoint.c
	$(CC) $(CFLAGS) -c ../MIDPOINT/midpoint.c -o midpoint.o

electrode.o				 : ../ELECTRODE/electrode.c
	$(CC) $(CFLAGS) -c ../ELECTRODE/electrode.c -o electrode.o

neighll.o					 : ../NEIGHLL/neighll.c
	$(CC) $(CFLAGS) -c ../NEIGHLL/neighll.c -o neighll.o

neighnode.o 			 : ../NEIGHNODE/neighnode.c
	$(CC) $(CFLAGS) -c ../NEIGHNODE/neighnode.c -o neighnode.o

node.o						 : ../NODE/node.c
	$(CC) $(CFLAGS) -c ../NODE/node.c -o node.o

omll.o 						 : ../OMLL/omll.c
	$(CC) $(CFLAGS) -c ../OMLL/omll.c -o omll.o

arbarray.o				 : ../ARBARRAY/arbarray.c
	$(CC) $(CFLAGS) -c ../ARBARRAY/arbarray.c -o arbarray.o

matrix.o 					 : ../MATRIX/matrix.c
	$(CC) $(CFLAGS) -c ../MATRIX/matrix.c -o matrix.o

matrix_linklist.o	 : ../MATRIX_LINKLIST/matrix_linklist.c
	$(CC) $(CFLAGS) -c ../MATRIX_LINKLIST/matrix_linklist.c -o matrix_linklist.o

clustersitenode.o  : ../CLUSTERSITENODE/clustersitenode.c
	$(CC) $(CFLAGS) -c ../CLUSTERSITENODE/clustersitenode.c -o clustersitenode.o

io.o 							 : ../IO/io.c
	$(CC) $(CFLAGS) -c ../IO/io.c -o io.o

.PHONY: clean
clean:
	$(RM) *.o
