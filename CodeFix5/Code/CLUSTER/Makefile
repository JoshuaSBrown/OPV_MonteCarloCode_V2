# Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=gnu99 -pedantic -g

# Targets
ALL : test_cluster test_cluster2

test_cluster : test_cluster.c  cluster.o node.o neighnode.o neighll.o electrode.o midpoint.o
	$(CC) $(CFLAGS) test_cluster.c cluster.o node.o neighnode.o neighll.o electrode.o midpoint.o -o test_cluster

test_cluster2 : test_cluster2.c  cluster.o node.o neighnode.o neighll.o electrode.o midpoint.o
	$(CC) $(CFLAGS) test_cluster2.c cluster.o node.o neighnode.o neighll.o electrode.o midpoint.o -o test_cluster2

cluster.o : cluster.c
	$(CC) $(CFLAGS) -c cluster.c -o cluster.o

node.o : ../NODE/node.c
	$(CC) $(CFLAGS) -c ../NODE/node.c -o node.o

neighll.o : ../NEIGHLL/neighll.c
	$(CC) $(CFLAGS) -c ../NEIGHLL/neighll.c -o neighll.o

neighnode.o : ../NEIGHNODE/neighnode.c
	$(CC) $(CFLAGS) -c ../NEIGHNODE/neighnode.c -o neighnode.o

electrode.o : ../ELECTRODE/electrode.c
	$(CC) $(CFLAGS) -c ../ELECTRODE/electrode.c -o electrode.o

midpoint.o : ../MIDPOINT/midpoint.c
	$(CC) $(CFLAGS) -c ../MIDPOINT/midpoint.c -o midpoint.o

.PHONY : clean
clean :
	$(RM) *.o *exe
