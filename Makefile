
CL = $(shell getconf LEVEL1_DCACHE_LINESIZE)

CFLAGS = -std=c++2a -Wall -Wextra -Wshadow -pedantic -march=native -DCACHE_LINE=$(CL)

.PHONY: clean

.DEFAULT: binoms

%/%.hpp:

binoms: binoms.cpp
	g++ $(CFLAGS) -DNDEBUG -Ofast -o binoms binoms.cpp

clean:
	rm -f binoms binomials build
