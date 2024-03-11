#include <iostream>
#include <bitset>

#include <sdsl/bit_vectors.hpp>

#include "h0_bv.hpp"


int main(int argc, char const *argv[]) {
    h0::internal::mults coder;
    if (argc > 1) {
        sdsl::bit_vector bv;
        bool opened = sdsl::load_from_file(bv, argv[1]);
        if (!opened) {
            std::cerr << "sdsl bv load failed" << std::endl;
            exit(1);
        } else {
            std::cout << "sdsl bv load done" << std::endl;
        }
        sdsl::rank_support_v rs(&bv);
        sdsl::select_support_mcl ss(&bv);
        h0::h0_bv h0bv(bv);

        if (argc > 2) {
            for (uint16_t i = 2; i < argc; ++i) {
                uint64_t idx = std::stoull(argv[i]);            
                std::cerr << "===================" << idx << "============================" << std::endl;
                uint64_t bv_val = bv.get_int((idx / 63) * 63, 63);
                std::cerr << std::bitset<64>(bv_val) << ", " << bv_val << " from " << ((idx / 63) * 63) << std::endl;
                auto enc = coder.encode(bv_val);
                std::cerr << " -> " << enc.first << ", " << enc.second << std::endl;
                uint64_t dec_val = coder.decode(enc.first, enc.second);
                std::cerr << std::bitset<64>(dec_val) << ", " << dec_val << std::endl;
                bool compval = h0bv.access(idx);
                std::cerr << "at " << idx << " (" << (idx % 63) << ") = " << ((bv_val >> (idx % 63)) & 1) 
                          << " or " << ((dec_val >> (idx % 63)) & 1) << " or " << compval << std::endl;
            }
            exit(0);
        }
        

        for (uint64_t i = 0; i < bv.size(); ++i) {
            bool cont = bv[i];
            bool res = h0bv.access(i);
            if (cont != res) {
                std::cerr << i << " should be " << cont << ", but is " << res << std::endl;
                exit(1);
            }
            auto c_count = rs.rank(i);
            auto r_count = h0bv.rank(i);
            if (c_count != r_count) {
                std::cerr << "rank(" << i << ") should be " << c_count << ", but is " << r_count << std::endl;
                exit(1);
            }
            if (i > 0 && (i % 10000000) == 0) {
                std::cout << i << "... ok\r" << std::flush;
            }
        }
        auto lim = rs(bv.size() - 1);
        for (uint64_t i = 1; i < lim; ++i) {
            auto cont = ss.select(i);
            auto res = h0bv.select(i);
            if (cont != res) {
                std::cerr << "select(" << i << ") should be " << cont << ", but is " << res << std::endl;
                exit(1);
            }
            if (i > 0 && (i % 10000000) == 0) {
                std::cout << i << "... ok              \r" << std::flush;
            }
        }
        std::cerr << "\nOK!?!?!?!?" << std::endl;
        exit(0);
    }

    uint64_t inp = 4287426570378804213;
    for (uint16_t i = 8; i < 64; i += 8) {
        uint64_t nv = inp & ((uint64_t(1) << i) - 1);
        auto r = coder.encode(nv);
        auto v = coder.decode(r.first, r.second);
        if (nv != v) {
            std::cerr << std::bitset<64>(nv) << "\n" << std::bitset<64>(v) << std::endl;
            std::cerr << nv << " -> " << r.first << ", " << r.second << " -> " << v << std::endl;
            exit(1);
        }    
    }
    auto r = coder.encode(inp);
    auto v = coder.decode(r.first, r.second);
    if (inp != v) {
        std::cerr << std::bitset<64>(inp) << "\n" << std::bitset<64>(v) << std::endl;
        std::cerr << inp << " -> " << r.first << ", " << r.second << " -> " << v << std::endl;
        exit(1);
    }
    inp = 9223367638808264703;
    r = coder.encode(inp);
    v = coder.decode(r.first, r.second);
    if (inp != v) {
        std::cerr << std::bitset<64>(inp) << "\n" << std::bitset<64>(v) << std::endl;
        std::cerr << inp << " -> " << r.first << ", " << r.second << " -> " << v << std::endl;
        exit(1);
    }

    inp = 720575940379279359;
    r = coder.encode(inp);
    v = coder.decode(r.first, r.second);
    if (inp != v) {
        std::cerr << std::bitset<64>(inp) << "\n" << std::bitset<64>(v) << std::endl;
        std::cerr << inp << " -> " << r.first << ", " << r.second << " -> " << v << std::endl;
        exit(1);
    }

    for (uint64_t i = 0; i < 100000000; ++i) {
        r = coder.encode(i);
        v = coder.decode(r.first, r.second);
        if (i != v) {
            std::cerr << i << " -> " << r.first << ", " << r.second << " -> " << v << std::endl;
            exit(1);
        }
    }
    inp = 0b101010;
    while (inp < (uint64_t(1) << 63)) {
        r = coder.encode(inp);
        v = coder.decode(r.first, r.second);
        if (inp != v) {
            std::cerr << inp << " -> " << r.first << ", " << r.second << " -> " << v << std::endl;
            exit(1);
        }
        inp <<= 1;
    }
    
}
