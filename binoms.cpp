#include <cstdint>
#include <iostream>
#include <bitset>
#include <random>

#include "h0_bv.hpp"

int main(int argc, char const *argv[]) {
    if (argc > 1) {
        uint16_t i = 0;
        uint64_t inp = std::stoull(argv[1]);
        std::cout << "\n" << std::bitset<64>(inp) << ", " << inp << std::endl;
        auto e = h0::internal::mults::encode(inp);
        std::cout << e.first << ", " << e.second << std::endl;
        auto dec = h0::internal::mults::decode(e.first, e.second);
        std::cout << std::bitset<64>(dec) << ", " << dec << std::endl;
        if (inp != dec) {
            uint64_t pci = __builtin_popcountll(inp);
            uint64_t pcd = __builtin_popcountll(dec);
            std::cerr << "a " << pci << "\nb " << pcd << "\n";
            uint64_t ppci = __builtin_popcountll(inp >> 32);
            uint64_t ppcd = __builtin_popcountll(dec >> 32);
            std::cerr << "a " << ppci << ", " << (pci - ppci) 
                      << "\nb " << ppcd << ", " << (pcd - ppcd) << "\na ";
            for (i = 3; i < 4; --i) {
                std::cerr << __builtin_popcountll((inp >> (i * 16)) & 0b1111111111111111);
                if (i) {
                    std::cerr << ", ";
                }
            }
            std::cerr << "\nb ";
            for (i = 3; i < 4; --i) {
                std::cerr << __builtin_popcountll((dec >> (i * 16)) & 0b1111111111111111);
                if (i) {
                    std::cerr << ", ";
                }
            }
            std::cerr << "\na ";
            for (i = 7; i < 8; --i) {
                std::cerr << __builtin_popcountll((inp >> (i * 8)) & 0b11111111);
                if (i) {
                    std::cerr << ", ";
                }
            }
            std::cerr << "\nb ";
            for (i = 7; i < 8; --i) {
                std::cerr << __builtin_popcountll((dec >> (i * 8)) & 0b11111111);
                if (i) {
                    std::cerr << ", ";
                }
            }
            std::cerr << std::endl;
            exit(1);
        }
        exit(0);
    }
    std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_int_distribution<uint64_t> dist(0);

    uint64_t i = 0;
    while (++i) {
        uint64_t inp = dist(mt);
        auto e = h0::internal::mults::encode(inp);
        auto dec = h0::internal::mults::decode(e.first, e.second);
        if (inp != dec) {
            std::cout << "\n" << std::bitset<64>(inp) << ", " << inp << std::endl;
            std::cout << std::bitset<64>(dec) << ", " << dec << std::endl; 
            exit(1);
        }
        if (i % 10000000 == 0) {
            std::cout << i << "\r" << std::flush;
        }
    }
    std::cout << std::endl;

    
    return 0;
}
