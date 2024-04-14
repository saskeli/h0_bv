#!/bin/bash

if [[ $# -lt 1 ]]; then
    echo "Illegal number of parameters" >&2
    exit 2
fi

lscpu | grep 'Model name'
uname -r

defs=(
    SDSL_15_32
    SDSL_24_32
    SDSL_31_32
    SDSL_63_32
    SDSL_256_32
    H0R_64_32
    H0I_63_32
    H0I_64_32
    SDSL_15_32_NOOPT
    SDSL_24_32_NOOPT
    SDSL_31_32_NOOPT
    SDSL_63_32_NOOPT
    H0R_64_32_NOOPT
    H0I_63_32_NOOPT
    H0I_64_32_NOOPT
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
    rnd-50.16MB
    WT-DNA-1GB
    WT-WEB-1GB
    bv-dump.bin
)

mkdir -p res

make bins

if [ $? -ne 0 ]
then
    echo "make failure... Terminating...."
    exit 1
fi

for DS in ${data[@]};
do
    echo "tests for ${DS}"
    for TYPE in ${defs[@]};
    do
        if [[ $# -ge 2 ]]; then
            if [[ ! $TYPE =~ $2 ]]; then
                continue
            fi
        fi
        echo "      ${TYPE}"
        ./${TYPE} ${1}${DS} 32 > res/${TYPE}_${DS}.txt
    done
done