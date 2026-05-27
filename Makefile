CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -Iinclude

.PHONY: all test clean

all: test

test: src/ab_testing.c tests/test_ab_testing.c
	@mkdir -p build
	$(CC) $(CFLAGS) src/ab_testing.c tests/test_ab_testing.c -lm -o build/test_ab_testing
	./build/test_ab_testing

clean:
	rm -rf build/
