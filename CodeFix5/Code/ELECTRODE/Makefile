# Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=gnu99 -pedantic -g

# Targets
ALL : test_electrode

test_electrode : test_electrode.c  electrode.o 
	$(CC) $(CFLAGS) test_electrode.c electrode.o -o test_electrode

electrode.o : electrode.c
	$(CC) $(CFLAGS) -c electrode.c -o electrode.o

.PHONY : clean
clean :
	$(RM) *.o *exe
