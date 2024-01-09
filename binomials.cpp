#include <iostream>
#include <iomanip>
#include <cstdint>

int main(int argc, char const *argv[]) {
    if (argc < 3) {
        std::cerr << "goan" << std::endl;
    }
    uint64_t n = std::stoull(argv[1]);
    uint64_t k = std::stoull(argv[2]);
    for (uint64_t i = 1; i <= k; ++i) {
        std::cout << std::setw(7) << i;
    }
    std::cout << std::endl;
    for (uint64_t nn = 1; nn <= n; ++nn) {
        uint64_t trg = std::min(k, nn);
        uint64_t current = nn;
        std::cout << std::setw(7) << current;
        uint64_t kk = 2; 
        while (kk <= trg) {
            current *= nn - kk + 1;
            current /= kk;
            std::cout << std::setw(7) << current;
            ++kk;
        }
        std::cout << std::endl;
    }
}
