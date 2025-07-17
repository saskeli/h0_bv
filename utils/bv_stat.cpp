#include <iostream>
#include <cstdint>
#include <cmath>
#include <array>
#include <vector>

#include "sdsl/bit_vectors.hpp"

constexpr uint64_t gcd(uint64_t a, uint64_t b) {
    if (a != b) {
        if (a < b) {
            return gcd(b - a, a);
        } else {
            return gcd(a - b, b);
        }
    }
    return a;
}

int main(int argc, char const *argv[]) {
    if (argc < 2) {
        std::cerr << "bv path is required!" << std::endl;
        exit(1);
    }
    sdsl::bit_vector bv;
    if (sdsl::load_from_file(bv, argv[1])) {
        sdsl::bit_vector::rank_1_type rs(&bv);
        uint64_t size = bv.size();
        uint64_t ones = rs.rank(size);
        std::cout << argv[1] << " loaded with " << size << " elements\n" << ones << " one bits " << double(ones) / size << std::endl;
        double p = ones;
        p /= size;
        double res = p * std::log2(1 / p) + (1 - p) * std::log2(1 / (1 - p));
        std::cout << "H0 entropy: " << res << std::endl;
        std::cout << "H0 bits: " << std::ceil(res * size) << std::endl;
        std::array a = {15, 24, 32, 64};
        for (uint64_t bs : a) {
            std::vector<uint64_t> bins = {1};
            for (uint64_t i = 1; i <= bs; ++i) {
                uint64_t g = gcd(bs - i + 1, i);
                uint64_t m = (bs - i + 1) / g;
                uint64_t ov = bins.back() / (i / g);
                ov *= m;
                bins.push_back(ov);
            }
            for (uint64_t i = 0; i < bins.size(); ++i) {
                bins[i] = std::ceil(std::log2(bins[i]));
            }
            uint64_t blocks = (bv.size() + bs) / bs;
            uint64_t uniform = 0;
            double mean_h0 = 0;
            uint64_t block_bits = 0;
            uint64_t payload_bits = 0;
            for (uint64_t trg = 0; trg < size; trg += bs) {
                uint64_t end = trg + bs > size ? size : trg + bs;
                ones = rs.rank(end) - rs.rank(trg);
                payload_bits += bins[ones];
                if (ones == 0 or ones == bs) {
                    ++uniform;
                } else {
                    p = ones;
                    p /= bs;
                    res = p * std::log2(1 / p) + (1 - p) * std::log2(1 / (1 - p));
                    mean_h0 += res;
                    block_bits += std::ceil(res * bs);
                }
            }
            mean_h0 /= blocks;
            double block_bits_per_block = block_bits;
            block_bits_per_block /= blocks;
            double bits_per_bit = double(block_bits) / size;
            uint64_t c_bits = blocks * std::ceil(std::log2(bs));
            std::cout << "Block size " << bs << ":\n" 
                      << "  " << blocks << " blocks\n"
                      << "  " << uniform << " " << double(uniform) / blocks << " uniform blocks\n"
                      << "  " << mean_h0 << " mean block h0 entropy\n"
                      << "  " << block_bits << " h0 block bits\n"
                      << "  " << block_bits_per_block << " mean h0 bits per block\n" 
                      << "  " << bits_per_bit << " h0 block bits/bit\n" 
                      << "  " << c_bits << " C bits\n" 
                      << "  " << payload_bits << " payload bits\n"
                      << "  " << double(c_bits + payload_bits) / size << " theoretical bits/bit" << std::endl;
        }
    }
    return 0;
}


