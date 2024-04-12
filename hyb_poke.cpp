#include <iostream>
#include <random>
#include <cstdint>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/hyb_vector.hpp>

int main(int argc, char const *argv[]) {
    std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_int_distribution<uint64_t> dist(0, ~uint64_t(0));

    sdsl::bit_vector bv(256);
    uint64_t counter = 0;

    while (true) {
        if (counter % 1000000 == 0) {
            std::cout << counter << "\r";
        }
        for (uint16_t i = 0; i < 4; ++i) {
            bv.set_int(i * 64, dist(mt), 64);
        }
        sdsl::hyb_vector<> hyb(bv);
        sdsl::hyb_vector::rank_1_type rs_hyb(&hyb);
        sdsl::bit_vector::rank_1_type rs_bv(&bv);
        for (uint16_t i = 0; i < 256; ++i) {
            if (bv[i] != hyb[i]) {
                for (uint16_t w = 0; w < 4; ++w) {
                    std::cerr << bv.get_int(w * 64, 64) << " ";
                }
                std::cerr << ", " << i << " access" << std::endl;
                exit(1);
            }
            if (rs_bv[i] != rs_hyb[i]) {
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