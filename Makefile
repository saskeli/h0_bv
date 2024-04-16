
CL = $(shell getconf LEVEL1_DCACHE_LINESIZE)

CFLAGS = -std=c++2b -Wall -Wextra -Wshadow -pedantic -march=native -DCACHE_LINE=$(CL)

PF = -Ofast -DNDEBUG

ISDSL = -isystem ~/include -L ~/lib -lsdsl

SDSL_INCLUDE = ~/include/sdsl/

.PHONY: clean bins

.DEFAULT: H0R_64_32_OP

%/%.hpp:

sdsl_comp: sdsl_comp.cpp h0_it.hpp internal.hpp h0_63.hpp h0_it.hpp h0_bv.hpp
	g++ $(CFLAGS) -DNDEBUG -Ofast -fconstexpr-ops-limit=1000000000 -isystem ~/include -L ~/lib -o sdsl_comp sdsl_comp.cpp -lsdsl

poke: internal.hpp poke.cpp h0_it.hpp h0_63.hpp h0_bv.hpp
	g++ $(CFLAGS) $(PF) -isystem ~/include -L ~/lib -o poke poke.cpp -lsdsl

hyb_poke: hyb_poke.cpp hyb_it.hpp hyb_256.hpp hyb_vanilla.hpp
	g++ $(CFLAGS) $(PF) -o hyb_poke hyb_poke.cpp $(ISDSL)

wdb_poke: h0_gap.hpp wdb_poke.cpp
	g++ $(CFLAGS) $(PF) -fconstexpr-ops-limit=1000000000 -isystem ~/include -L ~/lib -o wdb_poke wdb_poke.cpp -lsdsl

SDSL_15_32: rrr_time_and_space.cpp rrr_helper.hpp
	cp -f rrr_helper.hpp $(SDSL_INCLUDE)rrr_helper.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 -o SDSL_15_32 rrr_time_and_space.cpp $(ISDSL)

SDSL_24_32: rrr_time_and_space.cpp rrr_helper.hpp
	cp -f rrr_helper.hpp $(SDSL_INCLUDE)rrr_helper.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 -o SDSL_24_32 rrr_time_and_space.cpp $(ISDSL)

SDSL_31_32: rrr_time_and_space.cpp rrr_helper.hpp
	cp -f rrr_helper.hpp $(SDSL_INCLUDE)rrr_helper.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=31 -o SDSL_31_32 rrr_time_and_space.cpp $(ISDSL)

SDSL_63_32: rrr_time_and_space.cpp rrr_helper.hpp
	cp -f rrr_helper.hpp $(SDSL_INCLUDE)rrr_helper.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=63 -o SDSL_63_32 rrr_time_and_space.cpp $(ISDSL)

SDSL_64_32: rrr_time_and_space.cpp rrr_helper.hpp
	cp -f rrr_helper.hpp $(SDSL_INCLUDE)rrr_helper.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=64 -o SDSL_64_32 rrr_time_and_space.cpp $(ISDSL)

SDSL_256_32: rrr_time_and_space.cpp rrr_helper.hpp
	cp -f rrr_helper.hpp $(SDSL_INCLUDE)rrr_helper.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=256 -o SDSL_256_32 rrr_time_and_space.cpp $(ISDSL)

H0R_64_32: h0_bv.hpp internal.hpp  rrr_time_and_space.cpp 
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=64 -DHACK='"h0_bv.hpp"' -DCLASSNAME='h0::h0_bv<>' -o H0R_64_32 rrr_time_and_space.cpp $(ISDSL)

H0I_63_32: h0_63.hpp internal.hpp  rrr_time_and_space.cpp 
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=63 -DHACK='"h0_63.hpp"' -DCLASSNAME='h0::h0_bv<>' -o H0I_63_32 rrr_time_and_space.cpp $(ISDSL)

H0I_64_32: h0_it.hpp internal.hpp  rrr_time_and_space.cpp 
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=64 -DHACK='"h0_it.hpp"' -DCLASSNAME='h0::h0_bv<>' -o H0I_64_32 rrr_time_and_space.cpp $(ISDSL)

SDSL_15_32_NOOPT: rrr_time_and_space.cpp rrr_helper.hpp
	cp -f rrr_helper.hpp $(SDSL_INCLUDE)rrr_helper.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 -DRRR_NO_OPT -o SDSL_15_32_NOOPT rrr_time_and_space.cpp $(ISDSL)

SDSL_24_32_NOOPT: rrr_time_and_space.cpp rrr_helper.hpp
	cp -f rrr_helper.hpp $(SDSL_INCLUDE)rrr_helper.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 -DRRR_NO_OPT -o SDSL_24_32_NOOPT rrr_time_and_space.cpp $(ISDSL)

SDSL_31_32_NOOPT: rrr_time_and_space.cpp rrr_helper.hpp
	cp -f rrr_helper.hpp $(SDSL_INCLUDE)rrr_helper.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=31 -DRRR_NO_OPT -o SDSL_31_32_NOOPT rrr_time_and_space.cpp $(ISDSL)

SDSL_63_32_NOOPT: rrr_time_and_space.cpp rrr_helper.hpp
	cp -f rrr_helper.hpp $(SDSL_INCLUDE)rrr_helper.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=63 -DRRR_NO_OPT -o SDSL_63_32_NOOPT rrr_time_and_space.cpp $(ISDSL)

SDSL_64_32_NOOPT: rrr_time_and_space.cpp rrr_helper.hpp
	cp -f rrr_helper.hpp $(SDSL_INCLUDE)rrr_helper.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=64 -DRRR_NO_OPT -o SDSL_64_32_NOOPT rrr_time_and_space.cpp $(ISDSL)

H0R_64_32_NOOPT: h0_bv.hpp internal.hpp  rrr_time_and_space.cpp 
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=64 -DRRR_NO_OPT -DHACK='"h0_bv.hpp"' -DCLASSNAME='h0::h0_bv<>' -o H0R_64_32_NOOPT rrr_time_and_space.cpp $(ISDSL)

H0I_63_32_NOOPT: h0_63.hpp internal.hpp  rrr_time_and_space.cpp 
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=63 -DRRR_NO_OPT -DHACK='"h0_63.hpp"' -DCLASSNAME='h0::h0_bv<>' -o H0I_63_32_NOOPT rrr_time_and_space.cpp $(ISDSL)

H0I_64_32_NOOPT: h0_it.hpp internal.hpp  rrr_time_and_space.cpp 
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=64 -DRRR_NO_OPT -DHACK='"h0_it.hpp"' -DCLASSNAME='h0::h0_bv<>' -o H0I_64_32_NOOPT rrr_time_and_space.cpp $(ISDSL)

RRR_15_32: rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 -DRRR15 -o RRR_15_32 rrr_time_and_space.cpp $(ISDSL)

H0GAP_15_32_1: h0_gap.hpp internal.hpp  rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 -fconstexpr-ops-limit=1000000000 -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<15, 1>' -o H0GAP_15_32_1 rrr_time_and_space.cpp $(ISDSL)

H0GAP_15_32_15: h0_gap.hpp internal.hpp  rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 -fconstexpr-ops-limit=1000000000 -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<15>' -o H0GAP_15_32_15 rrr_time_and_space.cpp $(ISDSL)

H0GAP_15_32_24: h0_gap.hpp internal.hpp  rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 -fconstexpr-ops-limit=1000000000 -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<15, 24>' -o H0GAP_15_32_24 rrr_time_and_space.cpp $(ISDSL)

H0GAP_15_32_32: h0_gap.hpp internal.hpp  rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 -fconstexpr-ops-limit=1000000000 -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<15, 32>' -o H0GAP_15_32_32 rrr_time_and_space.cpp $(ISDSL)

H0GAP_15_32_64: h0_gap.hpp internal.hpp  rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 -fconstexpr-ops-limit=1000000000 -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<15, 64>' -o H0GAP_15_32_64 rrr_time_and_space.cpp $(ISDSL)

H0GAP_24_32_1: h0_gap.hpp internal.hpp  rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 -fconstexpr-ops-limit=5000000000 -fconstexpr-loop-limit=10000000 -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<24, 1>' -o H0GAP_24_32_1 rrr_time_and_space.cpp $(ISDSL)

H0GAP_24_32_15: h0_gap.hpp internal.hpp  rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 -fconstexpr-ops-limit=1000000000 -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<24, 15>' -o H0GAP_24_32_15 rrr_time_and_space.cpp $(ISDSL)

H0GAP_24_32_24: h0_gap.hpp internal.hpp  rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 -fconstexpr-ops-limit=1000000000 -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<24, 24>' -o H0GAP_24_32_24 rrr_time_and_space.cpp $(ISDSL)

H0GAP_24_32_32: h0_gap.hpp internal.hpp  rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 -fconstexpr-ops-limit=1000000000 -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<24, 32>' -o H0GAP_24_32_32 rrr_time_and_space.cpp $(ISDSL)

H0GAP_24_32_64: h0_gap.hpp internal.hpp  rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 -fconstexpr-ops-limit=1000000000 -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<24, 64>' -o H0GAP_24_32_64 rrr_time_and_space.cpp $(ISDSL)

H0WDBS_15_32: wdbs15.hpp internal.hpp  rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 -DHACK='"wdbs15.hpp"' -DCLASSNAME='h0::h0_wdb<>' -o H0WDBS_15_32 rrr_time_and_space.cpp $(ISDSL)

H0DBS_15_32: dbs15.hpp internal.hpp rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 -DHACK='"dbs15.hpp"' -DCLASSNAME='h0::h0_wdb<>' -o H0DBS_15_32 rrr_time_and_space.cpp $(ISDSL)

H0WDBS_24_32: wdbs24.hpp internal.hpp  rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 -DHACK='"wdbs24.hpp"' -DCLASSNAME='h0::h0_wdb<>' -o H0WDBS_24_32 rrr_time_and_space.cpp $(ISDSL)

H0DBS_24_32: dbs24.hpp internal.hpp rrr_time_and_space.cpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 -DHACK='"dbs24.hpp"' -DCLASSNAME='h0::h0_wdb<>' -o H0DBS_24_32 rrr_time_and_space.cpp $(ISDSL)

HYBSDSL_256_32: hyb_vanilla.hpp rrr_time_and_space.cpp
	cp -f hyb_vanilla.hpp $(SDSL_INCLUDE)hyb_vector.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=256 -DHYB -o HYBSDSL_256_32 rrr_time_and_space.cpp $(ISDSL)

HYBIT_256_32: hyb_it.hpp rrr_time_and_space.cpp
	cp -f hyb_it.hpp $(SDSL_INCLUDE)hyb_vector.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=256 -DHYB -o HYBIT_256_32 rrr_time_and_space.cpp $(ISDSL)

HYBRRR_256_32: hyb_256.hpp rrr_time_and_space.cpp
	cp -f hyb_256.hpp $(SDSL_INCLUDE)hyb_vector.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=256 -DHYB -o HYBRRR_256_32 rrr_time_and_space.cpp $(ISDSL)

bins: H0DBS_15_32 H0DBS_24_32 H0GAP_15_32_1 H0GAP_24_32_1 SDSL_256_32 SDSL_24_32_NOOPT SDSL_15_32_NOOPT SDSL_31_32_NOOPT SDSL_63_32_NOOPT SDSL_64_32_NOOPT H0R_64_32_NOOPT H0I_63_32_NOOPT H0I_64_32_NOOPT SDSL_24_32 SDSL_15_32 SDSL_31_32 SDSL_63_32 SDSL_64_32 H0R_64_32 H0I_63_32 H0I_64_32 RRR_15_32 H0GAP_15_32_15 H0GAP_15_32_24 H0GAP_15_32_32 H0GAP_15_32_64 H0GAP_24_32_15 H0GAP_24_32_24 H0GAP_24_32_32 H0GAP_24_32_64 H0WDBS_15_32 H0WDBS_24_32 HYBSDSL_256_32 HYBIT_256_32 HYBRRR_256_32

clean:
	cp -f hyb_vanilla.hpp $(SDSL_INCLUDE)hyb_vector.hpp
	rm -f poke wdb_poke sdsl_comp H0DBS_15_32 H0DBS_24_32 H0GAP_15_32_1 H0GAP_24_32_1 SDSL_256_32 SDSL_24_32_NOOPT SDSL_15_32_NOOPT SDSL_31_32_NOOPT SDSL_63_32_NOOPT SDSL_64_32_NOOPT H0R_64_32_NOOPT H0I_63_32_NOOPT H0I_64_32_NOOPT SDSL_24_32 SDSL_15_32 SDSL_31_32 SDSL_63_32 SDSL_64_32 H0R_64_32 H0I_63_32 H0I_64_32 RRR_15_32 H0GAP_15_32_15 H0GAP_15_32_24 H0GAP_15_32_32 H0GAP_15_32_64 H0GAP_24_32_15 H0GAP_24_32_24 H0GAP_24_32_32 H0GAP_24_32_64 H0WDBS_15_32 H0WDBS_24_32 HYBSDSL_256_32 HYBIT_256_32 HYBRRR_256_32
