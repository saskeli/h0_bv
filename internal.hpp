#include <immintrin.h>

#include <cstdint>
#include <array>
#include <vector>

namespace h0 {
namespace internal {

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

template <uint16_t n = 8>
constexpr std::array<uint64_t, n + 1> binoms() {
    std::array<uint64_t, n + 1> ret;
    ret[0] = 1;
    ret[1] = n;
    for (uint64_t k = 2; k <= n; ++k) {
        uint64_t g = gcd(n - k + 1, k);
        uint64_t m = (n - k + 1) / g;
        ret[k] = ret[k - 1] / (k / g);
        ret[k] *= m;
    }
    return ret;
}

constexpr uint32_t next(uint32_t v) {
    uint32_t off = __builtin_ctzl(v);
    uint32_t fix = __builtin_ctzl(~(v >> off)) - 1;
    v += uint32_t(1) << off;
    v |= (uint32_t(1) << fix) - 1;
    return v;
}

template <uint8_t k, uint8_t first>
constexpr std::array<uint8_t, k> byte_examples() {
    std::array<uint8_t, k> ret;
    ret[0] = first;
    for (uint8_t i = 1; i < k; ++i) {
        ret[i] = next(ret[i - 1]);
    }
    return ret;
}

template <std::array b0, std::array b1, std::array b2, std::array b3,
          std::array b4, std::array b5, std::array b6, std::array b7,
          std::array b8>
constexpr std::array<uint8_t, 256> byte_map() {
    std::array<uint8_t, 256> ret;
    for (uint16_t i = 0; i < b0.size(); ++i) {
        ret[b0[i]] = i;
    }
    for (uint16_t i = 0; i < b1.size(); ++i) {
        ret[b1[i]] = i;
    }
    for (uint16_t i = 0; i < b2.size(); ++i) {
        ret[b2[i]] = i;
    }
    for (uint16_t i = 0; i < b3.size(); ++i) {
        ret[b3[i]] = i;
    }
    for (uint16_t i = 0; i < b4.size(); ++i) {
        ret[b4[i]] = i;
    }
    for (uint16_t i = 0; i < b5.size(); ++i) {
        ret[b5[i]] = i;
    }
    for (uint16_t i = 0; i < b6.size(); ++i) {
        ret[b6[i]] = i;
    }
    for (uint16_t i = 0; i < b7.size(); ++i) {
        ret[b7[i]] = i;
    }
    for (uint16_t i = 0; i < b8.size(); ++i) {
        ret[b8[i]] = i;
    }
    return ret;
}

template <class A0, class A1, class A2, class A3, class A4, class A5, class A6,
          class A7, class A8>
constexpr auto concat(A0 a0, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7,
                      A8 a8) {
    std::array<uint8_t, 256> ret;
    size_t i = 0;
    std::copy_n(a0.begin(), a0.size(), ret.begin() + i);
    i += a0.size();
    std::copy_n(a1.begin(), a1.size(), ret.begin() + i);
    i += a1.size();
    std::copy_n(a2.begin(), a2.size(), ret.begin() + i);
    i += a2.size();
    std::copy_n(a3.begin(), a3.size(), ret.begin() + i);
    i += a3.size();
    std::copy_n(a4.begin(), a4.size(), ret.begin() + i);
    i += a4.size();
    std::copy_n(a5.begin(), a5.size(), ret.begin() + i);
    i += a5.size();
    std::copy_n(a6.begin(), a6.size(), ret.begin() + i);
    i += a6.size();
    std::copy_n(a7.begin(), a7.size(), ret.begin() + i);
    i += a7.size();
    std::copy_n(a8.begin(), a8.size(), ret.begin() + i);
    return ret;
}

template <class A0, class A1, class A2, class A3, class A4, class A5, class A6,
          class A7>
constexpr auto concat_count(A0 a0, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6,
                            A7 a7) {
    std::array<uint16_t, 9> ret;
    ret[0] = 0;
    ret[1] = ret[0] + a0.size();
    ret[2] = ret[1] + a1.size();
    ret[3] = ret[2] + a2.size();
    ret[4] = ret[3] + a3.size();
    ret[5] = ret[4] + a4.size();
    ret[6] = ret[5] + a5.size();
    ret[7] = ret[6] + a6.size();
    ret[8] = ret[7] + a7.size();
    return ret;
}

template <std::array<uint64_t, 9> byte_binoms = binoms()>
class byte_mapping {
   private:
    static const constexpr std::array<uint8_t, byte_binoms[0]> b0 =
        byte_examples<byte_binoms[0], 0b00000000>();
    static const constexpr std::array<uint8_t, byte_binoms[1]> b1 =
        byte_examples<byte_binoms[1], 0b00000001>();
    static const constexpr std::array<uint8_t, byte_binoms[2]> b2 =
        byte_examples<byte_binoms[2], 0b00000011>();
    static const constexpr std::array<uint8_t, byte_binoms[3]> b3 =
        byte_examples<byte_binoms[3], 0b00000111>();
    static const constexpr std::array<uint8_t, byte_binoms[4]> b4 =
        byte_examples<byte_binoms[4], 0b00001111>();
    static const constexpr std::array<uint8_t, byte_binoms[5]> b5 =
        byte_examples<byte_binoms[5], 0b00011111>();
    static const constexpr std::array<uint8_t, byte_binoms[6]> b6 =
        byte_examples<byte_binoms[6], 0b00111111>();
    static const constexpr std::array<uint8_t, byte_binoms[7]> b7 =
        byte_examples<byte_binoms[7], 0b01111111>();
    static const constexpr std::array<uint8_t, byte_binoms[8]> b8 =
        byte_examples<byte_binoms[8], 0b11111111>();
    static const constexpr std::array<uint8_t, 256> f_mapping =
        byte_map<b0, b1, b2, b3, b4, b5, b6, b7, b8>();
    static const constexpr std::array<uint8_t, 256> b_mapping =
        concat(b0, b1, b2, b3, b4, b5, b6, b7, b8);
    static const constexpr std::array<uint16_t, 9> bc_mapping =
        concat_count(b0, b1, b2, b3, b4, b5, b6, b7);

   public:
    static constexpr uint8_t f_byte(uint16_t k, uint64_t f) {
        return b_mapping[bc_mapping[k] + f];
    }

    static constexpr uint64_t get_f(uint8_t b) { return f_mapping[b]; }
};

template <class A, class B>
constexpr void fill(uint16_t t_idx, uint16_t s_idx, uint16_t hop, A& trg,
                    B& src) {
    trg[t_idx] = src[s_idx];
    if (hop == 0) {
        return;
    }
    fill(t_idx * 2, s_idx - hop, hop / 2, trg, src);
    fill(t_idx * 2 + 1, s_idx + hop, hop / 2, trg, src);
}

template <uint16_t b>
constexpr std::array<std::array<uint64_t, b / 2>, b + 1> f_lims() {
    std::array<std::array<uint64_t, b / 2 + 1>, b + 1> tmp;
    auto bins = binoms<b / 2>();
    for (uint16_t k = 0; k <= b; ++k) {
        tmp[k][0] = 0;
        for (uint16_t i = 0; i < b / 2; ++i) {
            if (k > b / 2 && i < (k - b / 2)) {
                tmp[k][i + 1] = 0;
            } else if (i > k) {
                tmp[k][i + 1] = tmp[k][i];
            } else {
                tmp[k][i + 1] = tmp[k][i];
                tmp[k][i + 1] += bins[i] * bins[k - i];
            }
        }
    }
    std::array<std::array<uint64_t, b / 2>, b + 1> ret;
    for (uint16_t k = 0; k <= b; ++k) {
        ret[k][0] = tmp[k][b / 2];
        fill(1, b / 4, b / 8, ret[k], tmp[k]);
    }
    return ret;
}

}  // namespace internal
}  //namespace h0