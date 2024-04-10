#!/bin/bash

lscpu | grep 'Model name'
uname -r

defs=(
    SDSL_15_32_OP
    SDSL_24_32_OP
    SDSL_31_32_OP
    SDSL_63_32_OP
    H0R_64_32_OP
    H0I_63_32_OP
    H0I_64_32_OP
    SDSL_15_32
    SDSL_24_32
    SDSL_31_32
    SDSL_63_32
    H0R_64_32
    H0I_63_32
    H0I_64_32
    RRR_15_32
    H0GAP_15_32_15
    H0GAP_15_32_24
    H0GAP_15_32_32
    H0GAP_15_32_64
    H0GAP_24_32_15
    H0GAP_24_32_24
    H0GAP_24_32_32
    H0GAP_24_32_64
    H0WDBS_15_32
    H0WDBS_24_32
    HYBSDSL_256_32 
    HYBIT_256_32 
    HYBRRR_256_32 
)

data=(
    rnd_50.16MB
    WT-DNA-1GB
    WT-WEB-1GB
)

mkdir -p res

make bins

for DS in ${data[@]};
do
    echo "tests for ${DS}"
    for TYPE in ${defs[@]};
    do
        echo "      ${TYPE}"
        ./${TYPE} ${1}${DS} 32 > res/${TYPE}_${DS}.txt
    done
done