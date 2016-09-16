# Variables
CC = gcc 
CFLAGS = -lm -Wall -Wextra -std=gnu99 -pedantic -g 

# Targets
ALL : test_functions

test_functions : test_functions.c functions.o matrix_linklist.o midpoint.o electrode.o neighnode.o neighll.o node.o montecarlo.o sitenode.o charge.o cluster.o matrix.o read.o 
	$(CC) $(CFLAGS) test_functions.c functions.o matrix_linklist.o midpoint.o electrode.o neighnode.o neighll.o node.o montecarlo.o sitenode.o charge.o cluster.o matrix.o read.o -o test_functions

functions.o : functions.c
	$(CC) $(CFLAGS) -c functions.c -o functions.o

montecarlo.o : ../MONTECARLO/montecarlo.c 
	$(CC) $(CFLAGS) -c ../MONTECARLO/montecarlo.c -o montecarlo.o

matrix.o : ../MATRIX/matrix.c
	$(CC) $(CFLAGS) -c ../MATRIX/matrix.c -o matrix.o

matrix_linklist.o : ../MATRIX_LINKLIST/matrix_linklist.c
	$(CC) $(CFLAGS) -c ../MATRIX_LINKLIST/matrix_linklist.c -o matrix_linklist.o

electrode.o : ../ELECTRODE/electrode.c
	$(CC) $(CFLAGS) -c ../ELECTRODE/electrode.c -o electrode.o

node.o : ../NODE/node.c
	$(CC) $(CFLAGS) -c ../NODE/node.c -o node.o

neighnode.o : ../NEIGHNODE/neighnode.c
	$(CC) $(CFLAGS) -c ../NEIGHNODE/neighnode.c -o neighnode.o

neighll.o : ../NEIGHLL/neighll.c
	$(CC) $(CFLAGS) -c ../NEIGHLL/neighll.c -o neighll.o

sitenode.o : ../SITENODE/sitenode.c
	$(CC) $(CFLAGS)  -c ../SITENODE/sitenode.c -o sitenode.o

cluster.o : ../CLUSTER/cluster.c
	$(CC) $(CFLAGS) -c ../CLUSTER/cluster.c -o cluster.o

midpoint.o : ../MIDPOINT/midpoint.c
	$(CC) $(CFLAGS) -c ../MIDPOINT/midpoint.c -o midpoint.o

charge.o : ../CHARGE/charge.c
	$(CC) $(CFLAGS) -c ../CHARGE/charge.c -o charge.o

read.o : ../PARAMETERS/read.c
	$(CC) $(CFLAGS) -c ../PARAMETERS/read.c -o read.o

.PHONY : clean
clean :
	$(RM) ../*/*.o
