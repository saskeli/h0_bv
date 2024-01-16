#include <iostream>
#include <sdsl/bit_vectors.hpp>

#include "h0_bv.hpp"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "input sdsl bv file required" << std::endl;
        exit(1);
    }
    sdsl::bit_vector bv;
    bool res = sdsl::load_from_file(bv, argv[1]);
    if (!res) {
        std::cerr << "sdsl bv load failed" << std::endl;
        exit(1);
    }
    sdsl::rank_support_v<> rs(&bv);
    sdsl::select_support_mcl<> ss(&bv);
    h0::h0_bv<> rrr_bv(bv);

    uint64_t tot = rs.rank(bv.size());

    std::cout << "created necessary structures with " << tot << " 1-bits\n"
              << "access..." << std::flush;
    for (uint64_t i = 0; i < bv.size(); ++i) {
        bool c = bv[i];
        bool r = rrr_bv.access(i);
        if (r != c) {
            std::cerr << "bv[" << i << "] = " << c << " != rrr_bv.access(" << i
                      << ") = " << r << std::endl;
            exit(1);
        }
    }
    std::cout << "ok\nrank..." << std::flush;
    for (uint64_t i = 0; i <= bv.size(); ++i) {
        uint64_t c = rs.rank(i);
        uint64_t r = rrr_bv.rank(i);
        if (r != c) {
            std::cerr << "rs.rank(" << i << ") = " << c << " != rrr_bv.rank("
                      << i << ") = " << r << std::endl;
            exit(1);
        }
    }
    std::cout << "ok\nselect..." << std::flush;
    for (uint64_t i = 1; i <= tot; ++i) {
        uint64_t c = ss.select(i);
        uint64_t r = rrr_bv.select(i);
        if (r != c) {
            std::cerr << "ss.select(" << i << ") = " << c
                      << " != rrr_bv.select(" << i << ") = " << r << std::endl;
            exit(1);
        }
    }
    std::cout << "ok" << std::endl;
    return 0;
}
