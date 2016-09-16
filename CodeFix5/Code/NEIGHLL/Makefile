# Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=gnu99 -pedantic -g

# Targets
ALL : test_neighll

test_neighll : test_neighll.c  neighll.o neighnode.o
	$(CC) $(CFLAGS) test_neighll.c neighll.o neighnode.o -o test_neighll

neighll.o : neighll.c
	$(CC) $(CFLAGS) -c neighll.c -o neighll.o

neighnode.o : ../NEIGHNODE/neighnode.c
	$(CC) $(CFLAGS) -c ../NEIGHNODE/neighnode.c -o neighnode.o

.PHONY : clean
clean :
	$(RM) *.o *exe
