#include <iostream>
#include <bitset>

#include "h0_gap.hpp"

int main() {
    h0::internal::mults<24, 24> decoder;
    std::cout << "made decoder" << std::endl;
    auto encoder = decoder.encode_array();
    std::cout << "made encodig array" << std::endl;

    uint32_t lim = uint32_t(1) << 24;

    /*for (uint32_t i = 0; i < 20; ++i) {
        std::cout << "8, " << i << " -> " << std::bitset<15>(decoder.decode(8, i)) << std::endl;
    }*/

    for (uint32_t i = 0; i < lim; ++i) {
        uint16_t C = __builtin_popcount(i);
        uint32_t off = encoder[i];
        auto v = decoder.decode<true>(C, off);
        if (v != i) {
            std::cerr << i << " != " << v << std::endl;
            std::cerr << std::bitset<31>(i) << "\n -> " << C << ", " << off << " -> \n" << std::bitset<31>(v) << std::endl;
            return 1;
        }
    }
    return 0;
}
