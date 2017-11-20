#Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -pedantic -g

#Targets
ALL : test_chargesitenode

test_chargesitenode : test_chargesitenode.c cluster.o electrode.o neighll.o midpoint.o neighnode.o node.o chargesitenode.o charge.o matrix.o matrix_linklist.o sitenode.o error.o 
	$(CC) $(CFLAGS) -o test_chargesitenode test_chargesitenode.c cluster.o electrode.o neighll.o midpoint.o neighnode.o node.o chargesitenode.o charge.o matrix.o matrix_linklist.o sitenode.o error.o

chargesitenode.o : chargesitenode.c
	$(CC) $(CFLAGS) -c chargesitenode.c

charge.o : ../CHARGE/charge.c
	$(CC) $(CFLAGS) -c ../CHARGE/charge.c

matrix.o : ../MATRIX/matrix.c
	$(CC) $(CFLAGS) -c ../MATRIX/matrix.c

matrix_linklist.o : ../MATRIX_LINKLIST/matrix_linklist.c
	$(CC) $(CFLAGS) -c ../MATRIX_LINKLIST/matrix_linklist.c

sitenode.o : ../SITENODE/sitenode.c
	$(CC) $(CFLAGS) -c ../SITENODE/sitenode.c

cluster.o : ../CLUSTER/cluster.c
	$(CC) $(CFLAGS) -c ../CLUSTER/cluster.c

neighnode.o : ../NEIGHNODE/neighnode.c
	$(CC) $(CFLAGS) -c ../NEIGHNODE/neighnode.c

node.o : ../NODE/node.c
	$(CC) $(CFLAGS) -c ../NODE/node.c

midpoint.o : ../MIDPOINT/midpoint.c
	$(CC) $(CFLAGS) -c ../MIDPOINT/midpoint.c

neighll.o : ../NEIGHLL/neighll.c
	$(CC) $(CFLAGS) -c ../NEIGHLL/neighll.c

electrode.o : ../ELECTRODE/electrode.c
	$(CC) $(CFLAGS) -c ../ELECTRODE/electrode.c

error.o : ../ERROR/error.c
	$(CC) $(CFLAGS) -c ../ERROR/error.c

.PHONY: clean
clean:
	$(RM) *.o

