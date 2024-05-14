
CL = $(shell getconf LEVEL1_DCACHE_LINESIZE)

CFLAGS = -std=c++2b -Wall -Wextra -Wshadow -pedantic -march=native -DCACHE_LINE=$(CL)

PF = -Ofast -DNDEBUG

ISDSL = -isystem ~/include -L ~/lib -lsdsl

INCLUDE = -I ./include

SDSL_INCLUDE = ~/include/sdsl/

RRR_TIMER = sdsl_hacks/rrr_time_and_space.cpp

CON_FLAGS = -fconstexpr-ops-limit=5000000000 -fconstexpr-loop-limit=10000000

BINARIES = build/H0GAP_24_32_7 build/H0GAP_15_32_7 build/H0LOO_15_32_1 build/H0LOO_15_32_15 \
build/H0LOO_24_32_1 build/H0LOO_24_32_24 build/H0DBS_15_32 build/H0DBS_24_32 build/H0GAP_15_32_1 \
build/H0GAP_24_32_1 build/SDSL_256_32 build/SDSL_24_32 build/SDSL_15_32 build/SDSL_31_32 \
build/SDSL_63_32 build/SDSL_64_32 build/H0R_64_32 build/H0I_63_32 build/H0I_64_32 build/RRR_15_32 \
build/H0GAP_15_32_15 build/H0GAP_15_32_24 build/H0GAP_15_32_32 build/H0GAP_15_32_64 \
build/H0GAP_24_32_15 build/H0GAP_24_32_24 build/H0GAP_24_32_32 build/H0GAP_24_32_64 \
build/H0WDBS_15_32 build/H0WDBS_24_32 build/HYBSDSL_256_32 build/HYBIT_256_32 build/HYBRRR_256_32

.PHONY: clean bins

.DEFAULT: bins

include/wdbs24.hpp: weightedDeBrujin.py
	python weightedDeBrujin.py 24 > include/wdbs24.hpp

include/wdbs15.hpp: weightedDeBrujin.py
	python weightedDeBrujin.py 15 > include/wdbs15.hpp

include/dbs24.hpp: weightedDeBrujin.py
	python weightedDeBrujin.py 24 full > include/dbs24.hpp

include/dbs15.hpp: weightedDeBrujin.py
	python weightedDeBrujin.py 15 full > include/dbs15.hpp

%/%.hpp:

build/SDSL_15_32: $(RRR_TIMER) sdsl_hacks/rrr_helper.hpp | build
	cp -f sdsl_hacks/rrr_helper.hpp $(SDSL_INCLUDE)rrr_helper.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 $(INCLUDE) -o build/SDSL_15_32 $(RRR_TIMER) $(ISDSL)

build/SDSL_24_32: $(RRR_TIMER) sdsl_hacks/rrr_helper.hpp | build
	cp -f sdsl_hacks/rrr_helper.hpp $(SDSL_INCLUDE)rrr_helper.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 $(INCLUDE) -o build/SDSL_24_32 $(RRR_TIMER) $(ISDSL)

build/SDSL_31_32: $(RRR_TIMER) sdsl_hacks/rrr_helper.hpp | build
	cp -f sdsl_hacks/rrr_helper.hpp $(SDSL_INCLUDE)rrr_helper.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=31 $(INCLUDE) -o build/SDSL_31_32 $(RRR_TIMER) $(ISDSL)

build/SDSL_63_32: $(RRR_TIMER) sdsl_hacks/rrr_helper.hpp | build
	cp -f sdsl_hacks/rrr_helper.hpp $(SDSL_INCLUDE)rrr_helper.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=63 $(INCLUDE) -o build/SDSL_63_32 $(RRR_TIMER) $(ISDSL)

build/SDSL_64_32: $(RRR_TIMER) sdsl_hacks/rrr_helper.hpp | build
	cp -f sdsl_hacks/rrr_helper.hpp $(SDSL_INCLUDE)rrr_helper.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=64 $(INCLUDE) -o build/SDSL_64_32 $(RRR_TIMER) $(ISDSL)

build/SDSL_256_32: $(RRR_TIMER) sdsl_hacks/rrr_helper.hpp | build
	cp -f sdsl_hacks/rrr_helper.hpp $(SDSL_INCLUDE)rrr_helper.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=256 $(INCLUDE) -o build/SDSL_256_32 $(RRR_TIMER) $(ISDSL)

build/H0R_64_32: include/h0_bv.hpp include/internal.hpp  $(RRR_TIMER)  | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=64 -DHACK='"h0_bv.hpp"' -DCLASSNAME='h0::h0_bv<>' $(INCLUDE) -o build/H0R_64_32 $(RRR_TIMER) $(ISDSL)

build/H0I_63_32: include/h0_63.hpp include/internal.hpp  $(RRR_TIMER)  | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=63 -DHACK='"h0_63.hpp"' -DCLASSNAME='h0::h0_bv<>' $(INCLUDE) -o build/H0I_63_32 $(RRR_TIMER) $(ISDSL)

build/H0I_64_32: include/h0_it.hpp include/internal.hpp  $(RRR_TIMER)  | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=64 -DHACK='"h0_it.hpp"' -DCLASSNAME='h0::h0_bv<>' $(INCLUDE) -o build/H0I_64_32 $(RRR_TIMER) $(ISDSL)

build/RRR_15_32: $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 -DRRR15 $(INCLUDE) -o build/RRR_15_32 $(RRR_TIMER) $(ISDSL)

build/H0GAP_15_32_1: include/h0_gap.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 $(CON_FLAGS) -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<15, 1>' $(INCLUDE) -o build/H0GAP_15_32_1 $(RRR_TIMER) $(ISDSL)

build/H0GAP_15_32_7: include/h0_gap.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 $(CON_FLAGS) -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<15, 7>' $(INCLUDE) -o build/H0GAP_15_32_7 $(RRR_TIMER) $(ISDSL)

build/H0GAP_15_32_15: include/h0_gap.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 $(CON_FLAGS) -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<15>' $(INCLUDE) -o build/H0GAP_15_32_15 $(RRR_TIMER) $(ISDSL)

build/H0GAP_15_32_24: include/h0_gap.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 $(CON_FLAGS) -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<15, 24>' $(INCLUDE) -o build/H0GAP_15_32_24 $(RRR_TIMER) $(ISDSL)

build/H0GAP_15_32_32: include/h0_gap.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 $(CON_FLAGS) -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<15, 32>' $(INCLUDE) -o build/H0GAP_15_32_32 $(RRR_TIMER) $(ISDSL)

build/H0GAP_15_32_64: include/h0_gap.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 $(CON_FLAGS) -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<15, 64>' $(INCLUDE) -o build/H0GAP_15_32_64 $(RRR_TIMER) $(ISDSL)

build/H0GAP_24_32_1: include/h0_gap.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 $(CON_FLAGS) -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<24, 1>' $(INCLUDE) -o build/H0GAP_24_32_1 $(RRR_TIMER) $(ISDSL)

build/H0GAP_24_32_7: include/h0_gap.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 $(CON_FLAGS) -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<24, 7>' $(INCLUDE) -o build/H0GAP_24_32_7 $(RRR_TIMER) $(ISDSL)

build/H0GAP_24_32_15: include/h0_gap.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 $(CON_FLAGS) -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<24, 15>' $(INCLUDE) -o build/H0GAP_24_32_15 $(RRR_TIMER) $(ISDSL)

build/H0GAP_24_32_24: include/h0_gap.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 $(CON_FLAGS) -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<24, 24>' $(INCLUDE) -o build/H0GAP_24_32_24 $(RRR_TIMER) $(ISDSL)

build/H0GAP_24_32_32: include/h0_gap.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 $(CON_FLAGS) -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<24, 32>' $(INCLUDE) -o build/H0GAP_24_32_32 $(RRR_TIMER) $(ISDSL)

build/H0GAP_24_32_64: include/h0_gap.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 $(CON_FLAGS) -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<24, 64>' $(INCLUDE) -o build/H0GAP_24_32_64 $(RRR_TIMER) $(ISDSL)

build/H0LOO_15_32_1: include/h0_gap.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 $(CON_FLAGS) -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<15, 1, 32, false>' $(INCLUDE) -o build/H0LOO_15_32_1 $(RRR_TIMER) $(ISDSL)
        
build/H0LOO_15_32_15: include/h0_gap.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 $(CON_FLAGS) -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<15, 15, 32, false>' $(INCLUDE) -o build/H0LOO_15_32_15 $(RRR_TIMER) $(ISDSL)

build/H0LOO_24_32_1: include/h0_gap.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 $(CON_FLAGS) -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<24, 1, 32, false>' $(INCLUDE) -o build/H0LOO_24_32_1 $(RRR_TIMER) $(ISDSL)

build/H0LOO_24_32_24: include/h0_gap.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 $(CON_FLAGS) -DHACK='"h0_gap.hpp"' -DCLASSNAME='h0::h0_gap<24, 24, 32, false>' $(INCLUDE) -o build/H0LOO_24_32_24 $(RRR_TIMER) $(ISDSL)

build/H0WDBS_15_32: include/wdbs15.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 -DHACK='"wdbs15.hpp"' -DCLASSNAME='h0::h0_wdb<>' $(INCLUDE) -o build/H0WDBS_15_32 $(RRR_TIMER) $(ISDSL)

build/H0DBS_15_32: include/dbs15.hpp include/internal.hpp $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=15 -DHACK='"dbs15.hpp"' -DCLASSNAME='h0::h0_wdb<>' $(INCLUDE) -o build/H0DBS_15_32 $(RRR_TIMER) $(ISDSL)

build/H0WDBS_24_32: include/wdbs24.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 -DHACK='"wdbs24.hpp"' -DCLASSNAME='h0::h0_wdb<>' $(INCLUDE) -o build/H0WDBS_24_32 $(RRR_TIMER) $(ISDSL)

build/H0DBS_24_32: include/dbs24.hpp include/internal.hpp $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=24 -DHACK='"dbs24.hpp"' -DCLASSNAME='h0::h0_wdb<>' $(INCLUDE) -o build/H0DBS_24_32 $(RRR_TIMER) $(ISDSL)

build/H0WDBS_31_32: include/wdbs31.hpp include/internal.hpp  $(RRR_TIMER) | build
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=31 -DHACK='"wdbs31.hpp"' -DCLASSNAME='h0::h0_wdb<>' $(INCLUDE) -o build/H0WDBS_31_32 $(RRR_TIMER) $(ISDSL)

build/HYBSDSL_256_32: sdsl_hacks/hyb_vanilla.hpp $(RRR_TIMER) | build
	cp -f sdsl_hacks/hyb_vanilla.hpp $(SDSL_INCLUDE)hyb_vector.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=256 -DHYB $(INCLUDE) -o build/HYBSDSL_256_32 $(RRR_TIMER) $(ISDSL)

build/HYBIT_256_32: sdsl_hacks/hyb_it.hpp $(RRR_TIMER) | build
	cp -f sdsl_hacks/hyb_it.hpp $(SDSL_INCLUDE)hyb_vector.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=256 -DHYB $(INCLUDE) -o build/HYBIT_256_32 $(RRR_TIMER) $(ISDSL)

build/HYBRRR_256_32: sdsl_hacks/hyb_256.hpp $(RRR_TIMER) | build
	cp -f sdsl_hacks/hyb_256.hpp $(SDSL_INCLUDE)hyb_vector.hpp
	g++ $(CFLAGS) $(PF) -DBLOCK_SIZE=256 -DHYB $(INCLUDE) -o build/HYBRRR_256_32 $(RRR_TIMER) $(ISDSL)

build: 
	mkdir build

bins: $(BINARIES)

clean:
	cp -f sdsl_hacks/hyb_vanilla.hpp $(SDSL_INCLUDE)hyb_vector.hpp
	rm -rf build
