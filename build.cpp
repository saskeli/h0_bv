#include <cstdint>
#include <bitset>
#include <iostream>
#include <chrono>
#include <array>
#include <algorithm>
#include <random>

#include "counters.hpp"

static const constexpr uint64_t MOD = 1000000007;

static constexpr uint64_t inv(uint64_t a) {
    return a <= 1 ? a : MOD - (MOD / a) * inv(MOD % a) % MOD;
}

template <uint16_t b>
static constexpr std::array<uint64_t, b + 1> factarr() {
    std::array<uint64_t, b + 1> arr;
    arr[0] = 1;
    for (uint32_t i = 1; i <= b; ++i) {
        arr[i] = arr[i - 1] * i;
        arr[i] %= MOD;
    }
    return arr;
}

template <uint16_t b>
static constexpr std::array<uint64_t, b + 1> invfactarr() {
    std::array<uint64_t, b + 1> arr = factarr<b>();
    for (uint32_t i = 0; i <= b; ++i) {
        arr[i] = inv(arr[i]);
    }
    return arr;
}

template <uint16_t b, std::array facts, std::array ifacts>
class BK_build {
   private:
    static_assert(b <= 32);
    uint64_t ncr(uint16_t n, uint16_t r) {
        uint64_t res = facts[n];
        res *= ifacts[r];
        res %= MOD;
        res *= ifacts[n - r];
        return res % MOD;
    }

    template <uint16_t sb>
    uint32_t tbuild(uint16_t k, uint64_t f, uint64_t olim) {
        if constexpr (sb == 1) {
            return k;
        } else {
            if (k == 1) {
                return k << f;
            }
            if (sb == k) {
                return (uint32_t(1) << k) - 1;
            }
            uint64_t r_sub = ncr(sb - 1, k - 1);
            uint64_t lim = olim - r_sub;
            const constexpr uint32_t m = uint32_t(1) << (sb - 1);
            uint32_t ind = f >= lim;
            r_sub = ind ? r_sub : lim;
            return (m * ind) | tbuild<sb - 1>(k - ind, f - lim * ind, r_sub);
        }
    }

   public:

    uint32_t build(uint16_t k, uint64_t f) {
        return tbuild<b>(k, f, ncr(b, k));
    }
};

template <class A>
void p(A& a) {
     
}

int main(int argc, char const *argv[]) {
    BK_build<32, factarr<32>(), invfactarr<32>()> bb;
    const constexpr uint32_t count = 601080390;
    uint32_t* ress = (uint32_t*)malloc(count * sizeof(uint32_t));
    for (uint32_t i = 0; i < count; ++i) {
        ress[i] = i;
    }
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(ress, ress+count, g);

    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::nanoseconds;

    Counters<1> counter;
    auto start = high_resolution_clock::now();
    for (uint32_t i = 0; i < count; ++i) {
        ress[i] = bb.build(8, ress[i]);
    }
    auto end = high_resolution_clock::now();
    auto a = counter.accumulate(0);

    std::cerr << double(duration_cast<nanoseconds>(end - start).count()) / count << std::endl;
    std::cerr << double(a[0]) / count << " cycles per query\n"
              << double(a[1]) / count << " retired instructions per query\n"
              << (double(a[1]) / a[0]) << " ipc\n"
              << double(a[2]) / count << " branch misspredictions per query\n"
              << double(a[3]) / count << " L1D misses per query" << std::endl;

    for (uint32_t i = 0; i < count; ++i) {
        std::cout << std::bitset<32>(ress[i]) << std::endl;
    }
    return 0;
}
