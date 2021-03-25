CC=gcc
CFLAGS=--std=c99 -Wall -pedantic -Wextra -DNDEBUG -O3

LIBS=

%: %.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
