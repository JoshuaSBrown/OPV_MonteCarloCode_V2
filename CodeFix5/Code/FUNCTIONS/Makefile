# Variables
CC = gcc 
CFLAGS = -lm -Wall -Wextra -std=gnu99 -pedantic -g 
INCLUDES += -I../SITENODE/ -I../CHARGE -I../MONTECARLO/


# Targets
ALL : test_functions

test_functions : test_functions.c functions.o montecarlo.o sitenode.o charge.o cluster.o matrix.o read.o 
	$(CC) $(CFLAGS) test_functions.c functions.o montecarlo.o sitenode.o charge.o cluster.o matrix.o read.o -o test_functions

functions.o : functions.c
	$(CC) $(CFLAGS) -c functions.c -o functions.o

montecarlo.o : ../CLUSTER/CLUSTERFUNCTIONS/MONTECARLO/montecarlo.c 
	$(CC) $(CFLAGS) -c ../CLUSTER/CLUSTERFUNCTIONS/MONTECARLO/montecarlo.c -o montecarlo.o

matrix.o : ../CLUSTER/CLUSTERFUNCTIONS/MATRIX/matrix.c
	$(CC) $(CFLAGS) -c ../CLUSTER/CLUSTERFUNCTIONS/MATRIX/matrix.c -o matrix.o

sitenode.o : ../CLUSTER/CLUSTERFUNCTIONS/SITENODE/sitenode.c
	$(CC) $(CFLAGS)  -c ../CLUSTER/CLUSTERFUNCTIONS/SITENODE/sitenode.c -o sitenode.o

cluster.o : ../CLUSTER/CLUSTERFUNCTIONS/DATASTRUCT/cluster.c
	$(CC) $(CFLAGS) -c ../CLUSTER/CLUSTERFUNCTIONS/DATASTRUCT/cluster.c -o cluster.o

charge.o : ../CHARGE/charge.c
	$(CC) $(CFLAGS) -c ../CHARGE/charge.c -o charge.o

read.o : ../PARAMETERS/read.c
	$(CC) $(CFLAGS) -c ../PARAMETERS/read.c -o read.o
#mem.o : ../MEM/mem.c
#	$(CC) $(CFLAGS) -c ../MEM/mem.c -o mem.o

.PHONY : clean
clean :
	$(RM) ../*/*.o
