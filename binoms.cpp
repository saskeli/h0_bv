
#include <array>
#include <bitset>
#include <cstdint>
#include <iostream>

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

   public:
    static constexpr uint8_t f_byte(uint16_t k, uint64_t f) {
        switch (k) {
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
        }
    }

    static constexpr uint64_t get_f(uint8_t b) { return f_mapping[b]; }
};

template <uint16_t n, std::array<uint64_t, n + 1> b64,
          std::array<uint64_t, (n / 2) + 1> b32,
          std::array<uint64_t, (n / 4) + 1> b16,
          std::array<uint64_t, (n / 8) + 1> b8>
class mults {
   private:
    template <uint16_t b>
    inline constexpr uint64_t encode(uint16_t k, uint64_t bv) {
        if (k == 0) {
            //std::cerr << b << " 0f: " << 0 << std::endl;
            return 0;
        } else if (k == 1) {
            //std::cerr << b << " 1f: " << __builtin_ctzll(bv) << std::endl;
            return __builtin_ctzll(bv);
        }
        if constexpr (b == 8) {
            //std::cerr << "byte f: " << byte_mapping<>::get_f(bv) << std::endl;
            return byte_mapping<>::get_f(bv);
        }
        const constexpr uint64_t mask = (uint64_t(1) << (b / 2)) - 1;
        uint64_t bvs = bv & mask;
        uint16_t ks = __builtin_popcountll(bvs);
        uint16_t kp = k - ks;
        uint64_t fp = encode<b / 2>(kp, bv >> (b / 2));
        uint64_t fs = encode<b / 2>(ks, bvs);
        uint64_t tot = 0;
        uint16_t start = k > (b / 2) ? k - (b / 2) : 0;
        for (uint16_t i = start; i < kp; ++i) {
            if constexpr (b == 64) {
                tot += b32[i] * b32[k - i];
            } else if constexpr (b == 32) {
                tot += b16[i] * b16[k - i];
            } else {
                //std::cerr << "tot " << tot << " -> ";
                tot += b8[i] * b8[k - i];
                //std::cerr << tot << std::endl;
            }
        }
        if constexpr (b == 64) {
            tot += fp * b32[ks];
        } else if constexpr (b == 32) {
            tot += fp * b16[ks];
        } else {
            //std::cerr << "tot " << tot << " inc to ";
            tot += fp * b8[ks];
            //std::cerr << tot << std::endl;
        }
        tot += fs;
        //std::cerr << b << " f: " << tot << std::endl;
        return tot;
    }

    template <uint16_t b>
    inline constexpr uint64_t decode(uint16_t k, uint64_t f) {
        if constexpr (b == 8) {
            return byte_mapping<>::f_byte(k, f);
        }
        uint16_t kp = k > (b / 2) ? k - (b / 2) : 0;
        uint64_t f_lim = 0;
        uint16_t end = k < (b / 2) ? k : b / 2;
        for (uint16_t i = kp; i < end; ++i) {
            uint64_t n_f_lim;
            if constexpr (b == 64) {
                n_f_lim = f_lim + b32[i] * b32[k - i];
            } else if constexpr (b == 32) {
                n_f_lim = f_lim + b16[i] * b16[k - i];
            } else {
                n_f_lim = f_lim + b8[i] * b8[k - i];
            }
            if (n_f_lim <= f) {
                ++kp;
                f_lim = n_f_lim;
            } else {
                break;
            }
        }
        uint16_t ks = k - kp;
        f -= f_lim;
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

   public:
    void print() {
        for (uint32_t i = 0; i < b64.size(); ++i) {
            std::cout << i << " " << b64[i] << std::endl;
        }
    }

    std::pair<uint16_t, uint64_t> encode(uint64_t bv) {
        /*std::cerr << "encoding " << bv << "(" << std::bitset<8>(bv >> 56) << " "
                  << std::bitset<8>(bv >> 48) << " " << std::bitset<8>(bv >> 40)
                  << " " << std::bitset<8>(bv >> 32) << " "
                  << std::bitset<8>(bv >> 24) << " " << std::bitset<8>(bv >> 16)
                  << " " << std::bitset<8>(bv >> 8) << " " << std::bitset<8>(bv)
                  << ")" << std::endl;*/

        uint16_t k = __builtin_popcountll(bv);
        return {k, encode<64>(k, bv)};
    }

    uint64_t decode(uint16_t k, uint64_t f) { return decode<64>(k, f); }
};

int main() {
    mults<64, binoms<64>(), binoms<32>(), binoms<16>(), binoms<8>()> m;
    uint64_t inp;
    uint16_t k;
    char c;
    /*for (uint64_t i = 1; i > 0; ++i) {
        auto r = m.encode(i);
        uint64_t v = m.decode(r.first, r.second);
        if (i != v) {
            std::cerr << i << " -> " << r.first << ", " << r.second << " -> " << v << std::endl;
            break;
        }
        if (i % 10000000 == 0) {
            std::cerr << i << std::endl;
        }
    }*/
    while (true) {
        std::cin >> c;
        if (c == 'e') {
            std::cin >> inp;
            auto r = m.encode(inp);
            std::cout << r.first << ", " << r.second << std::endl;
        } else if (c == 'd') {
            std::cin >> k >> inp;
            std::cout << m.decode(k, inp) << std::endl;
        } else {
            break;
        }
    }
    return 0;
}
