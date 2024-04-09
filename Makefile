
CL = $(shell getconf LEVEL1_DCACHE_LINESIZE)

CFLAGS = -std=c++2b -Wall -Wextra -Wshadow -pedantic -march=native -DCACHE_LINE=$(CL)

PF = -Ofast -DNDEBUG

ISDSL = -isystem ~/include -L ~/lib -lsdsl

.PHONY: clean bins

.DEFAULT: H0R_64_32_OP

%/%.hpp:

sdsl_comp: sdsl_comp.cpp h0_it.hpp internal.hpp h0_63.hpp h0_it.hpp h0_bv.hpp
	g++ $(CFLAGS) -DNDEBUG -Ofast -fconstexpr-ops-limit=1000000000 -isystem ~/include -L ~/lib -o sdsl_comp sdsl_comp.cpp -lsdsl

poke: internal.hpp poke.cpp h0_it.hpp h0_63.hpp h0_bv.hpp
	g++ $(CFLAGS) $(PF) -isystem ~/include -L ~/lib -o poke poke.cpp -lsdsl

wdb_poke: h0_gap.hpp wdb_poke.cpp
	g++ $(CFLAGS) $(PF) -fconstexpr-ops-limit=1000000000 -isystem ~/include -L ~/lib -o wdb_poke wdb_poke.cpp -lsdsl

SDSL_15_32_OP: rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 -o SDSL_15_32_OP rrr_time_and_space.cpp $(ISDSL)

SDSL_31_32_OP: rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=31 -o SDSL_31_32_OP rrr_time_and_space.cpp $(ISDSL)

SDSL_63_32_OP: rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=63 -o SDSL_63_32_OP rrr_time_and_space.cpp $(ISDSL)

H0R_64_32_OP: h0_bv.hpp rrr_time_and_space.cpp 
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=64 -DHACK='"h0_bv.hpp"' -DCLASSNAME='h0::h0_bv<>' -o H0R_64_32_OP rrr_time_and_space.cpp $(ISDSL)

H0I_63_32_OP: h0_63.hpp rrr_time_and_space.cpp 
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=63 -DHACK='"h0_63.hpp"' -DCLASSNAME='h0::h0_bv<>' -o H0I_63_32_OP rrr_time_and_space.cpp $(ISDSL)

H0I_64_32_OP: h0_it.hpp rrr_time_and_space.cpp 
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=64 -DHACK='"h0_it.hpp"' -DCLASSNAME='h0::h0_bv<>' -o H0I_64_32_OP rrr_time_and_space.cpp $(ISDSL)

SDSL_15_32: rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DRRR_NO_OPT -DBLOCK_SIZE=15 -o SDSL_15_32 rrr_time_and_space.cpp $(ISDSL)

SDSL_31_32: rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DRRR_NO_OPT -DBLOCK_SIZE=31 -o SDSL_31_32 rrr_time_and_space.cpp $(ISDSL)

SDSL_63_32: rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DRRR_NO_OPT -DBLOCK_SIZE=63 -o SDSL_63_32 rrr_time_and_space.cpp $(ISDSL)

H0R_64_32: h0_bv.hpp rrr_time_and_space.cpp 
	g++ $(CFLAGS) $(PF) -DRRR_NO_OPT -DBLOCK_SIZE=64 -DHACK='"h0_bv.hpp"' -DCLASSNAME='h0::h0_bv<>' -o H0R_64_32 rrr_time_and_space.cpp $(ISDSL)

H0I_63_32: h0_63.hpp rrr_time_and_space.cpp 
	g++ $(CFLAGS) $(PF) -DRRR_NO_OPT -DBLOCK_SIZE=63 -DHACK='"h0_63.hpp"' -DCLASSNAME='h0::h0_bv<>' -o H0I_63_32 rrr_time_and_space.cpp $(ISDSL)

H0I_64_32: h0_it.hpp rrr_time_and_space.cpp 
	g++ $(CFLAGS) $(PF) -DRRR_NO_OPT -DBLOCK_SIZE=64 -DHACK='"h0_it.hpp"' -DCLASSNAME='h0::h0_bv<>' -o H0I_64_32 rrr_time_and_space.cpp $(ISDSL)

RRR_15_32: rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 -DRRR15 -o RRR_15_32 rrr_time_and_space.cpp $(ISDSL)

H0GAP_15_32: h0_gap.hpp rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 -fconstexpr-ops-limit=1000000000 -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<15>' -o H0GAP_15_32 rrr_time_and_space.cpp $(ISDSL)

H0GAP_24_32: h0_gap.hpp rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 -fconstexpr-ops-limit=1000000000 -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<24>' -o H0GAP_24_32 rrr_time_and_space.cpp $(ISDSL)

H0WDBS_15_32: wdbs15.hpp rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 -DHACK='"wdbs15.hpp"' -DCLASSNAME='h0::h0_wdb<>' -o H0WDBS_15_32 rrr_time_and_space.cpp $(ISDSL)

H0WDBS_24_32: wdbs24.hpp rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 -DHACK='"wdbs24.hpp"' -DCLASSNAME='h0::h0_wdb<>' -o H0WDBS_24_32 rrr_time_and_space.cpp $(ISDSL)

bins: SDSL_15_32_OP SDSL_31_32_OP SDSL_63_32_OP H0R_64_32_OP H0I_63_32_OP H0I_64_32_OP SDSL_15_32 SDSL_31_32 SDSL_63_32 H0R_64_32 H0I_63_32 H0I_64_32 RRR_15_32 H0GAP_15_32 H0GAP_24_32 H0WDBS_15_32 H0WDBS_24_32

clean:
	rm -f SDSL_15_32_OP SDSL_31_32_OP SDSL_63_32_OP H0R_64_32_OP H0I_63_32_OP H0I_64_32_OP SDSL_15_32 SDSL_31_32 SDSL_63_32 H0R_64_32 H0I_63_32 H0I_64_32 RRR_15_32 H0GAP_15_32 H0GAP_24_32 H0WDBS_15_32 H0WDBS_24_32
