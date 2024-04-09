#!/bin/bash

defs=(
    "15 32 SDSL"
    "31 32 SDSL"
    "63 32 SDSL"
    "15 32 RRR15"
    "64 32 HACK=\"h0_bv.hpp\""
    "63 32 HACK=\"h0_63.hpp\""
    "64 32 HACK=\"h0_it.hpp\""
    "15 32 HACK=\"wdbs15.hpp\""
    "24 32 HACK=\"wdbs24.hpp\""
    "15 32 HACK=\"h0_gap.hpp\""
    "24 32 HACK=\"h0_gap.hpp\""
)

for 

g++ -std=c++2b -march=native -Ofast -isystem -DNDEBUG $defs ~/include -L ~/lib -o rrr_time_and_space rrr_time_and_space.cpp -lsdsl