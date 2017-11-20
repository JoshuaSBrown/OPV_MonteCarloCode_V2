# Variables
CC = gcc 
CFLAGS = -lm -Wall -Wextra -std=gnu99 -pedantic -g 

# Targets
ALL : test_io

test_io : test_io.c io.o  
	$(CC) $(CFLAGS) test_io.c io.o -o test_io

io.o : io.c
	$(CC) $(CFLAGS) -c io.c -o io.o

.PHONY : clean
clean :
	$(RM) *.o
