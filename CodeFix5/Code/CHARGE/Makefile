# Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -pedantic -g 

# Targets
ALL : test_charge

test_charge : test_charge.c charge.o 
	$(CC) $(CFLAGS) -o test_charge test_charge.c charge.o 

charge.o : charge.c
	$(CC) $(CFLAGS) -c charge.c

#mem.o : ../MEM/mem.c
#	$(CC) $(CFLAGS) -c ../MEM/mem.c

.PHONY: clean
clean:
	$(RM) *.o
