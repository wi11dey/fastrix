CC=gcc
CFLAGS=--std=c99 -Wall -pedantic -Wextra -DNDEBUG -O3

LIBS=

.PHONY: all
all: strassen

%: %.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
