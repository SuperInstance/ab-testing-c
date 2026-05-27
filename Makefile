CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -Iinclude
SRC = src/ab_testing.c
LIB = libabtesting.a

.PHONY: all lib test clean

all: lib

build:
	mkdir -p build

lib: $(LIB)

$(LIB): $(SRC) include/ab_testing.h | build
	$(CC) $(CFLAGS) -c $(SRC) -o build/ab_testing.o -lm
	ar rcs $@ build/ab_testing.o

test: test_runner
	./test_runner

test_runner: tests/test_ab_testing.c $(LIB) include/ab_testing.h
	$(CC) $(CFLAGS) tests/test_ab_testing.c -L. -labtesting -lm -o test_runner

clean:
	rm -f $(LIB) test_runner build/*.o
