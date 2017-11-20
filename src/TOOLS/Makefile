# Variables
CC = gcc
CFLAGS = -lm -Wall -Wextra -std=gnu99 -pedantic -g 

# Targets
ALL : run_analysis

run_analysis : run_analysis.c linklist2.o matrix.o analysis.o 
	$(CC) $(CFLAGS) run_analysis.c linklist2.o matrix.o analysis.o -o run_analysis

analysis.o : analysis.c linklist2.o matrix.o
	$(CC) $(CFLAGS) -c matrix.o linklist2.o analysis.c -o analysis.o

matrix.o : ../MATRIX/matrix.c
	$(CC) $(CFLAGS) -c ..//MATRIX/matrix.c -o matrix.o

linklist2.o : ../LINKLIST2/linklist2.c
	$(CC) $(CFLAGS) -c ../LINKLIST2/linklist2.c -o linklist2.o

.PHONY: clean
clean:
	$(RM) *.o
