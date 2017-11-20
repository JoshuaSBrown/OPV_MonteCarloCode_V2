# Variables
CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -pedantic -g
LDLIBS = -lm

# Targets
ALL : test_site

test_site : test_site.c sitenode.o 
	$(CC) $(CFLAGS) -o test_site test_site.c sitenode.o 

sitenode.o : sitenode.c 
	$(CC) $(CFLAGS) -c sitenode.c 

.PHONY: clean
clean:
	$(RM) *.o 

