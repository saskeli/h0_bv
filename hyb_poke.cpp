#include <iostream>
#include <random>
#include <cstdint>
#include <string>
#include <bitset>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/hyb_vector.hpp>

int main(int argc, char const *argv[]) {
    std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_int_distribution<uint64_t> dist(0, ~uint64_t(0));
    uint16_t n = 256;
    sdsl::bit_vector bv(n);

    uint64_t counter = 0;

    if (argc == 2) {
        bool opened = sdsl::load_from_file(bv, argv[1]);
        if (!opened) {
            std::cerr << "sdsl bv load failed" << std::endl;
            exit(1);
        } else {
            std::cout << "sdsl bv load done" << std::endl;
        }
        sdsl::hyb_vector<> hyb(bv);
        decltype(hyb)::rank_1_type rs_hyb(&hyb);
        decltype(bv)::rank_1_type rs_bv(&bv);
        for (uint64_t i = 0; i < bv.size(); ++i) {
            if (i % 10000000 == 0) {
                std::cout << i << "\r" << std::flush;
            }
            uint64_t cont = bv[i];
            uint64_t comp = hyb[i];
            if (cont != comp) {
                uint64_t idx = i / 256;
                uint64_t off = i % 256;
                std::cerr << "access(" << i << " -> (" << idx << ", " << off << ")) = " << comp << ", shoudl be " << cont << std::endl;
                for (uint16_t w = 0; w < 4; ++w) {
                    for (uint64_t woff = 0; woff < 64; ++woff) {
                        std::cerr << bv[idx * 256 + w * 64 + woff];
                    }
                    std::cerr << std::endl;
                }
                for (uint16_t w = 0; w < 4; ++w) {
                    std::cerr << bv.data()[idx * 4 + w] << " ";
                }
                std::cerr << std::endl;
                /*for (uint16_t w = 0; w < 4; ++w) {
                    for (uint64_t woff = 0; woff < 64; ++woff) {
                        std::cerr << hyb[idx * 256 + w * 64 + woff];
                    }
                    std::cerr << std::endl;
                }*/
                exit(1);
            }
            cont = rs_bv.rank(i);
            comp = rs_hyb.rank(i);
            if (cont != comp) {
                std::cerr << "rank(" << i << ") = " << comp << ", should be " << cont << std::endl;
                exit(1);
            }
        }
        std::cout << "\nOK!" << std::endl;
        exit(0);
    }

    if (argc > 1) {
        if (argc < 5) {
            std::cerr << "Go away!" << std::endl;
            exit(1);
        }
        for (uint16_t i = 0; i < 4; ++i) {
            uint64_t val = std::stoull(argv[i + 1]);
            std::cout << std::bitset<64>(val) << std::endl;
            for (uint16_t j = 0; j < n / 256; ++j) {
                bv.set_int((i + j * 4) * 64, val, 64);    
            }
        }
        sdsl::hyb_vector<> hyb(bv);
        decltype(hyb)::rank_1_type rs_hyb(&hyb);
        decltype(bv)::rank_1_type rs_bv(&bv);
        for (uint16_t i = 0; i < n; ++i) {
            uint16_t cont = bv[i];
            uint16_t comp = hyb[i];
            if (cont != comp) {
                std::cerr << "access(" << i << ") = " << comp << ", shoudl be " << cont << std::endl;
                exit(1);
            }
            cont = rs_bv.rank(i);
            comp = rs_hyb.rank(i);
            if (cont != comp) {
                std::cerr << "rank(" << i << ") = " << comp << ", should be " << cont << std::endl;
                exit(1);
            }
        }
        std::cout << "OK!" << std::endl;
        exit(0);
    }

    while (true) {
        if (counter % 10000 == 0) {
            std::cout << counter << "\r" << std::flush;
        }
        for (uint16_t i = 0; i < 4; ++i) {
            uint64_t v = dist(mt);
            for (uint16_t j = 0; j < 4; ++j) {
                bv.set_int((i + j * 4) * 64, v, 64);
            }
            
        }
        sdsl::hyb_vector<> hyb(bv);
        sdsl::hyb_vector<>::rank_1_type rs_hyb(&hyb);
        sdsl::bit_vector::rank_1_type rs_bv(&bv);
        for (uint16_t i = 0; i < n; ++i) {
            if (bv[i] != hyb[i]) {
                for (uint16_t w = 0; w < 4; ++w) {
                    std::cerr << bv.get_int(w * 64, 64) << " ";
                }
                std::cerr << ", " << i << " access" << std::endl;
                exit(1);
            }
            if (rs_bv.rank(i) != rs_hyb.rank(i)) {
                for (uint16_t w = 0; w < 4; ++w) {
                    std::cerr << bv.get_int(w * 64, 64) << " ";
                }
                std::cerr << ", " << i << " rank" << std::endl;
                exit(1);
            }
        }
        counter++;
    }
    return 0;
}