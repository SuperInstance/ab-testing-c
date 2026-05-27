CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2
AR = ar

SRC = src/ab_testing.c
OBJ = $(SRC:.c=.o)
LIB = libab_testing.a
HDR = include/ab_testing.h

.PHONY: all lib test clean

all: lib

lib: $(LIB)

$(LIB): $(OBJ)
	$(AR) rcs $@ $^

src/%.o: src/%.c $(HDR)
	$(CC) $(CFLAGS) -Iinclude -c -o $@ $<

test: $(LIB) tests/test_ab_testing.c $(HDR)
	$(CC) $(CFLAGS) -Iinclude -o test_ab_testing tests/test_ab_testing.c $(LIB) -lm
	./test_ab_testing

clean:
	rm -f src/*.o $(LIB) test_ab_testing
