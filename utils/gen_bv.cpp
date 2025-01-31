#include <random>
#include <string>

#include "sdsl/bit_vectors.hpp"

int main(int argc, char const *argv[]) {
    if (argc < 2) {
        std::cerr << "bv path is required!" << std::endl;
        exit(1);
    }
    uint64_t elems = 8;
    elems *= 1024;
    elems *= 1024;
    elems *= 1024;
    sdsl::bit_vector bv(elems);
    std::mt19937_64 gen(99999997);
    std::uniform_real_distribution<double> rand; 
    for (int k = -1; k >= -10; --k) {
        double cmp_val = std::pow(2.0, k);
        for (uint64_t i = 0; i < elems; ++i) {
            bv[i] = rand(gen) < cmp_val;
        }
        std::string w = argv[1];
        if (not w.ends_with('/')) {
            w.push_back('/');
        }
        w += "RND";
        w += std::to_string(k);
        w += ".bin";
        sdsl::store_to_file(bv, w);
    }
    return 0;
}
