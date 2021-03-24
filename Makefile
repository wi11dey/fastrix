CC=gcc
CFLAGS=--std=c99 -Wall -pedantic -Wextra -g

LIBS=

%: %.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
