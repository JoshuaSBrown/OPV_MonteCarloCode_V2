# Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=gnu99 -pedantic -g

# Targets
ALL : test_arbarray

test_arbarray : test_arbarray.c  arbarray.o midpoint.o omll.o cluster.o electrode.o neighnode.o neighll.o node.o error.o 
	$(CC) $(CFLAGS) test_arbarray.c arbarray.o midpoint.o omll.o cluster.o electrode.o neighnode.o neighll.o node.o error.o -o test_arbarray

arbarray.o : arbarray.c
	$(CC) $(CFLAGS) -c arbarray.c -o arbarray.o

midpoint.o : ../MIDPOINT/midpoint.c
	$(CC) $(CFLAGS) -c ../MIDPOINT/midpoint.c -o midpoint.o

omll.o : ../OMLL/omll.c
	$(CC) $(CFLAGS) -c ../OMLL/omll.c -o omll.o

cluster.o : ../CLUSTER/cluster.c
	$(CC) $(CFLAGS) -c ../CLUSTER/cluster.c -o cluster.o

neighnode.o : ../NEIGHNODE/neighnode.c
	$(CC) $(CFLAGS) -c ../NEIGHNODE/neighnode.c -o neighnode.o

neighll.o : ../NEIGHLL/neighll.c
	$(CC) $(CFLAGS) -c ../NEIGHLL/neighll.c -o neighll.o

node.o : ../NODE/node.c
	$(CC) $(CFLAGS) -c ../NODE/node.c -o node.o

electrode.o : ../ELECTRODE/electrode.c
	$(CC) $(CFLAGS) -c ../ELECTRODE/electrode.c -o electrode.o

error.o : ../ERROR/error.c
	$(CC) $(CFLAGS) -c ../ERROR/error.c -o error.o

.PHONY : clean
clean :
	$(RM) *.o *exe
