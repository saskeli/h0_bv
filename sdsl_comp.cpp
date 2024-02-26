#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include <chrono>
#include <random>
#include <vector>

#include "h0_it.hpp"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "input sdsl bv file required" << std::endl;
        exit(1);
    }
    uint64_t n = 100000000;
    sdsl::bit_vector bv;
    bool res = sdsl::load_from_file(bv, argv[1]);
    if (!res) {
        std::cerr << "sdsl bv load failed" << std::endl;
        exit(1);
    } else {
        std::cout << "sdsl bv load done" << std::endl;
    }
    sdsl::rank_support_v<> rs(&bv);
#ifndef TIME
    sdsl::select_support_mcl<> ss(&bv);
#else
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::nanoseconds;
#endif
    h0::h0_bv<> rrr_bv(bv);
    uint64_t tot = rs.rank(bv.size());
    std::cout << "created necessary " << bv.size() << "-bit structures with " << tot << " 1-bits" << std::endl;

    std::vector<uint64_t> q;
    std::random_device rd;  
    std::mt19937 gen(argc > 2 ? std::stoull(argv[2]) : rd()); 
    std::uniform_int_distribution<uint64_t> distrib(0, bv.size() - 1);
 
    for (uint64_t i = 0; i < n; ++i) {
        q.push_back(distrib(gen));
    }

    std::cout << "access..." << std::flush;
#ifdef TIME
    uint64_t checksum = 0;
    auto s = high_resolution_clock::now();
#endif
    for (uint64_t i = 0; i < n; ++i) {
#ifndef TIME
        bool c = bv[q[i]];
#endif
        bool r = rrr_bv.access(q[i]);
#ifndef TIME
        if (r != c) {
            std::cerr << "sdsl.access[" << q[i] << "] = " << c << " != h0_bv.access(" << q[i]
                      << ") = " << r << std::endl;
            exit(1);
        }
#else
        checksum += r;
#endif
    }
#ifdef TIME
    auto e = high_resolution_clock::now();
    std::cout << double(duration_cast<nanoseconds>(e - s).count()) / n << "ns\t" << checksum << "\nrank..." << std::flush;
    checksum = 0;
    s = high_resolution_clock::now();
#else
    std::cout << "ok\nrank..." << std::flush;
#endif
    for (uint64_t i = 0; i < n; ++i) {
#ifndef TIME
        uint64_t c = rs.rank(q[i]);
#endif
        uint64_t r = rrr_bv.rank(q[i]);
#ifndef TIME
        if (r != c) {
            std::cerr << "sdsl.rank(" << q[i] << ") = " << c << " != h0_bv.rank("
                      << q[i] << ") = " << r << std::endl;
            exit(1);
        }
#else 
        checksum += r;
#endif
    }
#ifdef TIME
    e = high_resolution_clock::now();
    std::cout << double(duration_cast<nanoseconds>(e - s).count()) / n << "ns\t" << checksum << std::endl;
#else
    std::cout << "ok" << std::endl;
#endif
    
    distrib = std::uniform_int_distribution<uint64_t>(1, tot);
 
    for (uint64_t i = 0; i < n; ++i) {
        q[i] = distrib(gen);
    }
    std::cout << "select..." << std::flush;
#ifdef TIME
    checksum = 0;
    s = high_resolution_clock::now();
#endif
    for (uint64_t i = 1; i < n; ++i) {
#ifndef TIME
        uint64_t c = ss.select(q[i]);
#endif
        uint64_t r = rrr_bv.select(q[i]);
#ifndef TIME
        if (r != c) {
            std::cerr << "sdsl.select(" << q[i] << ") = " << c
                      << " != h0_bv.select(" << q[i] << ") = " << r << std::endl;
            exit(1);
        }
#else 
        checksum += r;
#endif
    }
#ifdef TIME
    e = high_resolution_clock::now();
    std::cout << double(duration_cast<nanoseconds>(e - s).count()) / n << "ns\t" << checksum << std::endl;
#else
    std::cout << "ok" << std::endl;
#endif
    return 0;
}
