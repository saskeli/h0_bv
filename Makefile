
CL = $(shell getconf LEVEL1_DCACHE_LINESIZE)

CFLAGS = -std=c++2b -Wall -Wextra -Wshadow -pedantic -march=native -DCACHE_LINE=$(CL)

.PHONY: clean

.DEFAULT: binoms

%/%.hpp:

binoms: binoms.cpp h0_bv.hpp
	g++ $(CFLAGS) -DNDEBUG -Ofast -o binoms binoms.cpp

sdsl_comp: sdsl_comp.cpp h0_it.hpp internal.hpp
	g++ $(CFLAGS) -DNDEBUG -Ofast -isystem ~/include -L ~/lib -o sdsl_comp sdsl_comp.cpp -lsdsl

poke: internal.hpp poke.cpp h0_it.hpp
	g++ $(CFLAGS) -Ofast -DNDEBUG -isystem ~/include -L ~/lib -o poke poke.cpp -lsdsl

clean:
	rm -f binoms binomials build
