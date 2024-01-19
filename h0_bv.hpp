#include <immintrin.h>

#include <algorithm>
#include <array>
#include <bitset>
#include <cstdint>
#include <iostream>
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
        /*switch (k) {
            case 0:
                return b0[f];
            case 1:
                return b1[f];
            case 2:
                return b2[f];
            case 3:
                return b3[f];
            case 4:
                return b4[f];
            case 5:
                return b5[f];
            case 6:
                return b6[f];
            case 7:
                return b7[f];
            default:
                return b8[f];
        }*/
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

class mults {
   private:
    static const constexpr uint16_t n = 64;
    static const constexpr std::array<uint64_t, n / 2 + 1> b32 = binoms<32>();
    static const constexpr std::array<uint64_t, n / 4 + 1> b16 = binoms<16>();
    static const constexpr std::array<uint64_t, n / 8 + 1> b8 = binoms<8>();
    static const constexpr std::array<std::array<uint64_t, n / 2>, n + 1>
        f_lim64 = f_lims<64>();
    static const constexpr std::array<std::array<uint64_t, n / 4>, n / 2 + 1>
        f_lim32 = f_lims<32>();
    static const constexpr std::array<std::array<uint64_t, n / 8>, n / 4 + 1>
        f_lim16 = f_lims<16>();

    template <uint16_t b>
    inline static constexpr uint64_t encode(uint16_t k, uint64_t bv) {
        if constexpr (b == 8) {
            // std::cerr << "byte f: " << byte_mapping<>::get_f(bv) <<
            // std::endl;
            return byte_mapping<>::get_f(bv);
        }
        if (k == 0) {
            // std::cerr << b << " 0f: " << 0 << std::endl;
            return 0;
        } else if (k == 1) {
            // std::cerr << b << " 1f: " << __builtin_ctzll(bv) << std::endl;
            return __builtin_ctzll(bv);
        } else if (k == b) {
            return 0;
        } else if (k == b - 1) {
            return b - __builtin_ctzll(~bv) - 1;
        }
        const constexpr uint64_t mask = (uint64_t(1) << (b / 2)) - 1;
        uint64_t bvs = bv & mask;
        uint16_t ks = __builtin_popcountll(bvs);
        uint16_t kp = k - ks;
        uint64_t fp = encode<b / 2>(kp, bv >> (b / 2));
        uint64_t fs = encode<b / 2>(ks, bvs);
        uint64_t tot = 0;
        uint16_t start = k > (b / 2) ? k - (b / 2) : 0;
        for (uint16_t i = start; i < b; ++i) {
            if (i == kp) {
                break;
            }
            if constexpr (b == 64) {
                tot += b32[i] * b32[k - i];
            } else if constexpr (b == 32) {
                tot += b16[i] * b16[k - i];
            } else {
                // std::cerr << "tot " << tot << " -> ";
                tot += b8[i] * b8[k - i];
                // std::cerr << tot << std::endl;
            }
        }

        if constexpr (b == 64) {
            tot += fp * b32[ks];
        } else if constexpr (b == 32) {
            tot += fp * b16[ks];
        } else {
            // std::cerr << "tot " << tot << " inc to ";
            tot += fp * b8[ks];
            // std::cerr << tot << std::endl;
        }
        tot += fs;
        // std::cerr << b << " f: " << tot << std::endl;
        return tot;
    }

    template <uint16_t b>
    inline static constexpr uint16_t kp_f(uint16_t k, uint64_t& f) {
        uint16_t kp = 0;
        uint16_t i = 1;
        uint64_t f_lim = 0;
        if constexpr (b == 64) {
            uint64_t fl_v = f_lim64[k][i];
            bool r = fl_v <= f;
            f_lim = r ? fl_v : f_lim;
            i = i * 2 + r;
            kp += r * b / 4;
            fl_v = f_lim64[k][i];
            r = fl_v <= f;
            f_lim = r ? fl_v : f_lim;
            i = i * 2 + r;
            kp += r * b / 8;
            fl_v = f_lim64[k][i];
            r = fl_v <= f;
            f_lim = r ? fl_v : f_lim;
            i = i * 2 + r;
            kp += r * b / 16;
            fl_v = f_lim64[k][i];
            r = fl_v <= f;
            f_lim = r ? fl_v : f_lim;
            i = i * 2 + r;
            kp += r * b / 32;
            fl_v = f_lim64[k][i];
            r = fl_v <= f;
            f_lim = r ? fl_v : f_lim;
            kp += r;
            fl_v = f_lim64[k][0];
            r = f >= fl_v;
            kp += r;
            f_lim = r ? fl_v : f_lim;
        } else if constexpr (b == 32) {
            uint64_t fl_v = f_lim32[k][i];
            bool r = fl_v <= f;
            i = i * 2 + r;
            kp += r * b / 4;
            f_lim = r ? fl_v : f_lim;
            fl_v = f_lim32[k][i];
            r = fl_v <= f;
            i = i * 2 + r;
            kp += r * b / 8;
            f_lim = r ? fl_v : f_lim;
            fl_v = f_lim32[k][i];
            r = fl_v <= f;
            i = i * 2 + r;
            kp += r * b / 16;
            f_lim = r ? fl_v : f_lim;
            fl_v = f_lim32[k][i];
            r = fl_v <= f;
            kp += r;
            f_lim = r ? fl_v : f_lim;
            fl_v = f_lim32[k][0];
            r = f >= fl_v;
            kp += r;
            f_lim = r ? fl_v : f_lim;
        } else {
            uint64_t fl_v = f_lim16[k][i];
            bool r = fl_v <= f;
            i = i * 2 + r;
            kp += r * b / 4;
            f_lim = r ? fl_v : f_lim;
            fl_v = f_lim16[k][i];
            r = fl_v <= f;
            i = i * 2 + r;
            kp += r * b / 8;
            f_lim = r ? fl_v : f_lim;
            fl_v = f_lim16[k][i];
            r = fl_v <= f;
            kp += r;
            f_lim = r ? fl_v : f_lim;
            fl_v = f_lim16[k][0];
            r = f >= fl_v;
            kp += r;
            f_lim = r ? fl_v : f_lim;
        }
        f -= f_lim;
        return kp;
    }

    template <uint16_t b>
    inline static constexpr uint64_t decode(uint16_t k, uint64_t f) {
        if constexpr (b == 8) {
            return byte_mapping<>::f_byte(k, f);
        }
        if (k == b) {
            if constexpr (b == 64) {
                return ~uint64_t(0);
            } else {
                return (uint64_t(1) << b) - 1;
            }
        } else if (k == 0) {
            return 0;
        } else if (k == 1) {
            return uint64_t(1) << f;
        } else if (k == b - 1) {
            if (b == 64) {
                return (~uint64_t(0)) ^ (uint64_t(1) << (b - f - 1));
            } else {
                return ((uint64_t(1) << b) - 1) ^ (uint64_t(1) << (b - f - 1));
            }
        }
        uint16_t kp = kp_f<b>(k, f);
        uint16_t ks = k - kp;
        uint64_t md;
        if constexpr (b == 64) {
            md = b32[ks];
        } else if constexpr (b == 32) {
            md = b16[ks];
        } else {
            md = b8[ks];
        }
        uint64_t ret = decode<b / 2>(kp, f / md) << (b / 2);
        ret |= decode<b / 2>(ks, f % md);
        return ret;
    }

    template <uint16_t b>
    inline static constexpr uint64_t rank(uint16_t k, uint64_t f, uint16_t i) {
        if constexpr (b == 8) {
            auto by = byte_mapping<>::f_byte(k, f);
            return __builtin_popcount(by & ((uint16_t(1) << i) - 1));
        }
        if (i >= b) {
            return k;
        } else if (k == 0) {
            return 0;
        } else if (k == b) {
            return i;
        } else if (k == 1) {
            return f < i;
        } else if (k == b - 1) {
            return i - uint64_t((b - f - 1) < i);
        }
        uint16_t kp = kp_f<b>(k, f);
        uint16_t ks = k - kp;
        uint64_t md;
        if constexpr (b == 64) {
            md = b32[ks];
        } else if constexpr (b == 32) {
            md = b16[ks];
        } else {
            md = b8[ks];
        }
        if (i == b / 2) {
            return ks;
        } else if (i > b / 2) {
            return ks + rank<b / 2>(kp, f / md, i - b / 2);
        }
        return rank<b / 2>(ks, f % md, i);
    }

    template <uint16_t b>
    inline static constexpr bool access(uint16_t k, uint64_t f, uint16_t i) {
        if constexpr (b == 8) {
            auto by = byte_mapping<>::f_byte(k, f);
            return (by >> i) & uint8_t(1);
        }
        if (k == 0) {
            return 0;
        } else if (k == b) {
            return 1;
        } else if (k == 1) {
            return f == i;
        } else if (k == b - 1) {
            return (b - f - 1) != i;
        }
        uint16_t kp = kp_f<b>(k, f);
        uint16_t ks = k - kp;
        uint64_t md;
        if constexpr (b == 64) {
            md = b32[ks];
        } else if constexpr (b == 32) {
            md = b16[ks];
        } else {
            md = b8[ks];
        }
        if (i < b / 2) {
            return access<b / 2>(ks, f % md, i);
        }
        return access<b / 2>(kp, f / md, i - b / 2);
    }

    template <uint16_t b>
    inline static constexpr uint64_t select(uint16_t k, uint64_t f,
                                            uint16_t s) {
        if constexpr (b == 8) {
            uint32_t by = byte_mapping<>::f_byte(k, f);
            uint16_t pos = uint16_t(1) << (s - 1);
            return __builtin_ctz(_pdep_u32(pos, by));
        }
        if (k == 0) {
            return 0;
        } else if (k == b) {
            return s - 1;
        } else if (k == 1) {
            return f;
        } else if (k == b - 1) {
            return s - 1 + ((b - f - 1) < s);
        }
        uint16_t kp = kp_f<b>(k, f);
        uint16_t ks = k - kp;
        uint64_t md;
        if constexpr (b == 64) {
            md = b32[ks];
        } else if constexpr (b == 32) {
            md = b16[ks];
        } else {
            md = b8[ks];
        }
        if (ks < s) {
            return b / 2 + select<b / 2>(kp, f / md, s - ks);
        }
        return select<b / 2>(ks, f % md, s);
    }

   public:
    inline static constexpr std::pair<uint16_t, uint64_t> encode(uint64_t bv) {
        /*std::cerr << "encoding " << bv << "(" << std::bitset<8>(bv >> 56) << "
           "
                  << std::bitset<8>(bv >> 48) << " " << std::bitset<8>(bv >> 40)
                  << " " << std::bitset<8>(bv >> 32) << " "
                  << std::bitset<8>(bv >> 24) << " " << std::bitset<8>(bv >> 16)
                  << " " << std::bitset<8>(bv >> 8) << " " << std::bitset<8>(bv)
                  << ")" << std::endl;*/

        uint16_t k = __builtin_popcountll(bv);
        return {k, encode<64>(k, bv)};
    }

    inline static constexpr uint64_t decode(uint16_t k, uint64_t f) {
        return decode<64>(k, f);
    }

    inline static constexpr uint64_t rank(uint16_t k, uint64_t f, uint16_t i) {
        return rank<64>(k, f, i);
    }

    inline static constexpr bool access(uint16_t k, uint64_t f, uint16_t i) {
        // return access<64>(k, f, i);
        // std::cerr << "access(" << k << ", " << f << ", " << i << ")" <<
        // std::endl;
        if (k == 0) [[unlikely]] {
            return 0;
        }
        if (k == n) [[unlikely]] {
            return 1;
        }
        if (k == 1) [[unlikely]] {
            return f == i;
        }
        if (k == n - 1) [[unlikely]] {
            return (63 - f) != i;
        }
        uint16_t kp = kp_f<n>(k, f);
        uint16_t ks = k - kp;
        uint64_t md = b32[ks];
        bool r = i < n / 2;
        bool l = !r;
        k = ks * r + kp * l;
        f = (f % md) * r + (f / md) * l;
        i -= l * (n / 2);
        // std::cerr << "      (" << k << ", " << f << ", " << i << ")" <<
        // std::endl;
        /*if (k == 0) [[unlikely]] {
            return 0;
        }
        if (k == n / 2) [[unlikely]] {
            return 1;
        }
        if (k == 1) [[unlikely]] {
            return f == i;
        }
        if (k == n / 2 - 1) [[unlikely]] {
            return (31 - f) != i;
        }*/
        kp = kp_f<n / 2>(k, f);
        ks = k - kp;
        md = b16[ks];
        r = i < n / 4;
        l = !r;
        k = ks * r + kp * l;
        f = (f % md) * r + (f / md) * l;
        i -= l * (n / 4);
        // std::cerr << "      (" << k << ", " << f << ", " << i << ")" <<
        // std::endl;
        /*if (k == 0) [[unlikely]] {
            return 0;
        }
        if (k == n / 4) [[unlikely]] {
            return 1;
        }
        if (k == 1) [[unlikely]] {
            return f == i;
        }
        if (k == n / 4 - 1) [[unlikely]] {
            return (15 - f) != i;
        }*/
        kp = kp_f<n / 4>(k, f);
        ks = k - kp;
        md = b8[ks];
        r = i < n / 8;
        l = !r;
        k = ks * r + kp * l;
        f = (f % md) * r + (f / md) * l;
        i -= l * (n / 8);
        // std::cerr << "      (" << k << ", " << f << ", " << i << ")" <<
        // std::endl;
        auto by = byte_mapping<>::f_byte(k, f);
        return (by >> i) & uint8_t(1);
    }

    inline static constexpr uint64_t select(uint16_t k, uint64_t f,
                                            uint16_t s) {
        return select<64>(k, f, s);
    }
};

class v_ints {
   private:
    std::vector<uint64_t> data;

   public:
    void append(uint64_t v, uint64_t off, uint16_t w) {
        uint64_t mod = off % 64;
        uint64_t div = off / 64;
        if (div + 1 >= data.size()) {
            data.push_back(0);
        }
        data[div] |= v << mod;
        if (mod + w > 64) {
            data[div + 1] = v >> (64 - mod);
        }
    }

    uint64_t read(uint64_t off, uint16_t w) const {
        uint64_t mod = off % 64;
        uint64_t div = off / 64;
        uint64_t val = data[div] >> mod;
        if (mod + w > 64) {
            val |= data[div + 1] << (64 - mod);
        }
        val &= (uint64_t(1) << w) - 1;
        return val;
    }

    uint64_t size() const { return data.size(); }

    uint64_t bytes_size() const { return sizeof(v_ints) + size() * 8; }
};

template <uint16_t b>
const constexpr std::array<uint64_t, b + 1> f_widhts() {
    std::array<uint64_t, b + 1> ret;
    auto bins = binoms<b>();
    for (uint16_t i = 0; i <= b; ++i) {
        ret[i] = 64 - __builtin_clzll(bins[i] - 1);
    }
    return ret;
}

}  // namespace internal

template <uint16_t k = 32>
class h0_bv {
   private:
    static const constexpr internal::mults coder = internal::mults();
    static const constexpr auto widths = internal::f_widhts<64>();
    static const constexpr uint16_t block_width = 64;
    static const constexpr uint16_t k_width = 7;
    std::vector<std::pair<uint64_t, uint64_t>> meta;
    internal::v_ints data;
    uint64_t size_;
    uint64_t sum_;

   public:
    template <class bv_type>
    h0_bv(const bv_type& bv) : meta(), data(), sum_(0) {
        const uint64_t* bv_data = bv.data();
        size_ = bv.size();
        uint64_t block = 0;
        uint64_t offset = 0;
        uint64_t blocks = size_ / block_width + (size_ % block_width ? 1 : 0);
        while (blocks--) {
            if (block % k == 0) {
                meta.push_back({offset, sum_});
            }
            uint64_t elem = bv_data[block++];
            auto enc = coder.encode(elem);
            // std::cerr << std::bitset<64>(elem) << " -> " << enc.first << ", "
            // << enc.second << std::endl;
            sum_ += enc.first;
            data.append(enc.first, offset, k_width);
            offset += k_width;
            data.append(enc.second, offset, widths[enc.first]);
            // std::cerr << "                -> " << widths[enc.first] <<
            // std::endl;
            offset += widths[enc.first];
        }
    }

    bool access(uint64_t i) const {
        uint64_t block = i / block_width;
        uint64_t s_block = block / k;
        uint16_t b_offset = block % k;
        i %= block_width;
        uint64_t offset = meta[s_block].first;
        for (uint64_t j = 0; j < k; ++j) {
            if (b_offset == j) {
                break;
            }
            auto bk = data.read(offset, k_width);
            offset += k_width + widths[bk];
        }
        auto bk = data.read(offset, k_width);
        offset += k_width;
        auto f = data.read(offset, widths[bk]);
        return coder.access(bk, f, i);
    }

    uint64_t rank(uint64_t i) const {
        uint64_t block = i / block_width;
        uint64_t s_block = block / k;
        uint16_t b_offset = block % k;
        i %= block_width;
        uint64_t offset = meta[s_block].first;
        uint64_t prefix_rank = meta[s_block].second;
        for (uint64_t j = 0; j < k; ++j) {
            if (b_offset == j) {
                break;
            }
            auto bk = data.read(offset, k_width);
            prefix_rank += bk;
            offset += k_width + widths[bk];
        }
        auto bk = data.read(offset, k_width);
        offset += k_width;
        auto f = data.read(offset, widths[bk]);
        return prefix_rank + coder.rank(bk, f, i);
    }

    uint64_t select(uint64_t s) const {
        uint64_t a = 0;
        uint64_t b = meta.size() - 1;
        while (a < b) {
            uint64_t m = (a + b + 1) / 2;
            if (meta[m].second >= s) {
                b = m - 1;
            } else {
                a = m;
            }
        }
        uint64_t offset = meta[a].first;
        s -= meta[a].second;
        uint64_t res = k * a * block_width;
        uint16_t bk = 0;
        uint64_t f = 0;
        for (uint16_t i = 0; i < k; ++i) {
            bk = data.read(offset, k_width);
            offset += k_width;
            if (bk >= s) {
                f = data.read(offset, widths[bk]);
                break;
            }
            offset += widths[bk];
            res += block_width;
            s -= bk;
        }
        return res + coder.select(bk, f, s);
    }

    uint64_t size() const { return size_; }

    uint64_t bytes_size() const {
        return sizeof(h0_bv) + meta.size() * 16 + data.bytes_size();
    }
};

}  // namespace h0