CC       ?= gcc
CXX      ?= g++
CFLAGS   ?= -O2 -std=c99 -Wall -Wextra -pedantic -Werror
CXXFLAGS ?= -O2 -std=c++11 -Wall -Wextra -pedantic -Werror
LDLIBS    = -lm

all: test

test: test_c test_c_alt test_cpp
	./test_c
	./test_c_alt
	./test_cpp

test_c: test.c angle16.h
	$(CC) $(CFLAGS) -o $@ test.c $(LDLIBS)

# Same tests with the non-default implementation selected for the public API.
test_c_alt: test.c angle16.h
	$(CC) $(CFLAGS) -DANGLE16_TO_FLOAT_BITS=0 -DANGLE16_FROM_FLOAT_BITS=1 -o $@ test.c $(LDLIBS)

test_cpp: test_cpp.cpp angle16.h
	$(CXX) $(CXXFLAGS) -o $@ test_cpp.cpp $(LDLIBS)

bench: bench.c angle16.h
	$(CC) $(CFLAGS) -o $@ bench.c $(LDLIBS)

clean:
	rm -f test_c test_c_alt test_cpp bench

.PHONY: all test clean
