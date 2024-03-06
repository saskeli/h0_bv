#include <immintrin.h>

#include <algorithm>
#include <array>
#include <bitset>
#include <cstdint>
#include <iostream>
#include <vector>

#include <sdsl/int_vector.hpp>
#include "internal.hpp"

namespace h0 {
namespace internal {
class mults {
   private:
    static const constexpr uint16_t n = 64;
    static const constexpr std::array<uint64_t, n + 1> b64 = binoms<64>();
    static const constexpr std::array<uint64_t, n - 8 + 1> b56 = binoms<56>();
    static const constexpr std::array<uint64_t, n - 16 + 1> b48 = binoms<48>();
    static const constexpr std::array<uint64_t, n - 24 + 1> b40 = binoms<40>();
    static const constexpr std::array<uint64_t, n - 32 + 1> b32 = binoms<32>();
    static const constexpr std::array<uint64_t, n - 40 + 1> b24 = binoms<24>();
    static const constexpr std::array<uint64_t, n - 48 + 1> b16 = binoms<16>();
    static const constexpr std::array<uint64_t, n - 56 + 1> b8 = binoms<8>();
    static const constexpr std::array<std::array<uint64_t, 9>, n + 1>
        f_lim64 = f_lims_it<64>();
    static const constexpr std::array<std::array<uint64_t, 9>, n - 8 + 1>
        f_lim56 = f_lims_it<56>();
    static const constexpr std::array<std::array<uint64_t, 9>, n - 16 + 1>
        f_lim48 = f_lims_it<48>();
    static const constexpr std::array<std::array<uint64_t, 9>, n - 24 + 1>
        f_lim40 = f_lims_it<40>();
    static const constexpr std::array<std::array<uint64_t, 9>, n - 32 + 1>
        f_lim32 = f_lims_it<32>();
    static const constexpr std::array<std::array<uint64_t, 9>, n - 40 + 1>
        f_lim24 = f_lims_it<24>();
    static const constexpr std::array<std::array<uint64_t, 9>, n - 48 + 1>
        f_lim16 = f_lims_it<16>();

    template <uint16_t b>
    inline static constexpr uint64_t encode(uint16_t k, uint64_t bv) {
        if constexpr (b <= 8) {
            // std::cerr << "byte f: " << byte_mapping<>::get_f(bv) <<
            // std::endl;
            return byte_mapping<>::get_f(bv);
        } else {
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
            const constexpr uint64_t mask = (uint64_t(1) << (b - 8)) - 1;
            uint64_t bvs = bv & mask;
            uint16_t ks = __builtin_popcountll(bvs);
            uint16_t kp = k - ks;
            uint64_t fp = byte_mapping<>::get_f(bv >> (b - 8));
            uint64_t fs = encode<b - 8>(ks, bvs);
            uint64_t tot = 0;
            if constexpr (b == 64) {
                tot += f_lim64[k][kp];
                tot += fp * b56[ks];
            } else if constexpr (b == 56) {
                tot += f_lim56[k][kp];
                tot += fp * b48[ks];
            } else if constexpr (b == 48) {
                tot += f_lim48[k][kp];
                tot += fp * b40[ks];
            } else if constexpr (b == 40) {
                tot += f_lim40[k][kp];
                tot += fp * b32[ks];
            } else if constexpr (b == 32) {
                tot += f_lim32[k][kp];
                tot += fp * b24[ks];
            } else if constexpr (b == 24) {
                tot += f_lim24[k][kp];
                tot += fp * b16[ks];
            } else {
                tot += f_lim16[k][kp];
                tot += fp * b8[ks];
            }

            tot += fs;
            // std::cerr << b << " f: " << tot << std::endl;
            return tot;
        }
    }

    template<class T>
    inline static constexpr void bs(T& fl, uint16_t& kp, uint64_t& f_lim, uint64_t& f) {
        kp = 0;
        for (uint16_t i = 1; i <= 8; ++i) {
            kp += f >= fl[i];
        }
        f_lim = fl[kp];
        f -= f_lim;
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
        uint16_t kp;
        uint64_t f_lim;

        if (f == 0) {
            if (k == n) {
                return ~uint64_t(0);
            }
            return (uint64_t(1) << k) - 1;
        }
        if (f == b64[k] - 1) {
            return ~((uint64_t(1) << (n - k)) - 1);
        }
        if (k == 1) {
            return uint64_t(1) << f;
        }
        if (k == (n - 1)) {
            return ~(uint64_t(1) << (63 - f));
        }

        bs(f_lim64[k], kp, f_lim, f);
        k -= kp;
        uint64_t md = b56[k];
        uint64_t res = byte_mapping<>::f_byte(kp, f / md);
        f %= md;

        bs(f_lim56[k], kp, f_lim, f);
        k -= kp;
        md = b48[k];
        res = (res << 8) | byte_mapping<>::f_byte(kp, f / md);
        f %= md;

        bs(f_lim48[k], kp, f_lim, f);
        k -= kp;
        md = b40[k];
        res = (res << 8) | byte_mapping<>::f_byte(kp, f / md);
        f %= md;

        bs(f_lim40[k], kp, f_lim, f);
        k -= kp;
        md = b32[k];
        res = (res << 8) | byte_mapping<>::f_byte(kp, f / md);
        f %= md;

        bs(f_lim32[k], kp, f_lim, f);
        k -= kp;
        md = b24[k];
        res = (res << 8) | byte_mapping<>::f_byte(kp, f / md);
        f %= md;

        bs(f_lim24[k], kp, f_lim, f);
        k -= kp;
        md = b16[k];
        res = (res << 8) | byte_mapping<>::f_byte(kp, f / md);
        f %= md;

        bs(f_lim16[k], kp, f_lim, f);
        k -= kp;
        md = b8[k];
        res = (res << 8) | byte_mapping<>::f_byte(kp, f / md);

        res = (res << 8) | byte_mapping<>::f_byte(k, f % md);
        return res;
    }

    inline static constexpr uint64_t rank(uint16_t k, uint64_t f, uint16_t i) {
        return __builtin_popcountll(decode(k, f) & ((uint64_t(1) << i) - 1));
    }

    inline static constexpr bool access(uint16_t k, uint64_t f, uint16_t i) {
        return (decode(k, f) >> i) & uint64_t(1);
    }

    inline static uint64_t select(uint16_t k, uint64_t f,
                                            uint16_t s) {
        return 63 - __builtin_clzll(_pdep_u64(uint64_t(1) << (s - 1), decode(k, f)));
    }
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
} // namespace internal

template <uint16_t k = 32>
class h0_bv {
   private:
    static const constexpr internal::mults coder = internal::mults();
    static const constexpr auto widths = internal::f_widhts<64>();
    static const constexpr uint16_t block_width = 64;
    static const constexpr uint16_t k_width = 7;
    sdsl::int_vector<> meta_point;
    sdsl::int_vector<> meta_rank;
    sdsl::int_vector<> data_typ;
    sdsl::bit_vector data_val;
    uint64_t size_;
    uint64_t sum_;

   public:
    h0_bv(const sdsl::bit_vector& bv) : sum_(0) {
        const uint64_t* bv_data = bv.data();
        size_ = bv.size();
        uint64_t blocks = size_ / block_width + (size_ % block_width ? 1 : 0);
        uint64_t block = 0;
        uint64_t offset = 0;
        while (blocks--) {
            uint64_t elem = bv_data[block++];
            uint16_t p = __builtin_popcountll(elem);
            sum_ += p;
            offset += widths[p];
        }

        uint64_t v_siz = size_ / (k * block_width);
        v_siz += (size_ % (k * block_width)) > 0;
        data_typ = sdsl::int_vector<>(size_ / block_width + (size_ % block_width > 0), 0, k_width);
        data_val = sdsl::bit_vector(offset);
        meta_point = sdsl::int_vector<>(v_siz, 0, 64 - __builtin_clzll(size_));
        meta_rank = sdsl::int_vector<>(v_siz, 0, 64 - __builtin_clzll(sum_));

        block = 0;
        sum_ = 0;
        offset = 0;
        blocks = size_ / block_width + (size_ % block_width ? 1 : 0);
        uint64_t s_block = 0;
        while (blocks--) {
            if (block % k == 0) {
                assert(meta_point.size() > s_block);
                meta_point[s_block] = offset;
                assert(meta_rank.size() > s_block);
                meta_rank[s_block++] = sum_;
            }
            uint64_t elem = bv_data[block];
            auto enc = coder.encode(elem);
            // std::cerr << std::bitset<64>(elem) << " -> " << enc.first << ", "
            // << enc.second << std::endl;
            sum_ += enc.first;
            assert(data_typ.size() > block);
            data_typ[block++] = enc.first;
            data_val.set_int(offset, enc.second, widths[enc.first]);
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
        uint64_t offset = meta_point[s_block];
        for (uint64_t j = 0; j < k; ++j) {
            if (b_offset == j) {
                break;
            }
            auto bk = data_typ[s_block * k + j];
            offset += widths[bk];
        }
        auto bk = data_typ[block];
        auto f = data_val.get_int(offset, widths[bk]);
        return coder.access(bk, f, i);
    }

    uint64_t rank(uint64_t i) const {
        uint64_t block = i / block_width;
        uint64_t s_block = block / k;
        uint16_t b_offset = block % k;
        i %= block_width;
        uint64_t offset = meta_point[s_block];
        uint64_t prefix_rank = meta_rank[s_block];
        for (uint64_t j = 0; j < k; ++j) {
            if (b_offset == j) {
                break;
            }
            auto bk = data_typ[s_block * k + j];
            prefix_rank += bk;
            offset += widths[bk];
        }
        auto bk = data_typ[block];
        auto f = data_val.get_int(offset, widths[bk]);
        return prefix_rank + coder.rank(bk, f, i);
    }

    uint64_t select(uint64_t s) const {
        uint64_t a = 0;
        uint64_t b = meta_point.size() - 1;
        while (a < b) {
            uint64_t m = (a + b + 1) / 2;
            if (meta_rank[m] >= s) {
                b = m - 1;
            } else {
                a = m;
            }
        }
        uint64_t offset = meta_point[a];
        s -= meta_rank[a];
        uint64_t res = k * a * block_width;
        uint16_t bk = 0;
        uint64_t f = 0;
        uint64_t block = a * k;
        for (uint16_t i = 0; i < k; ++i) {
            bk = data_typ[block++];
            if (bk >= s) {
                f = data_val.get_int(offset, widths[bk]);
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
        return sizeof(h0_bv) + sdsl::size_in_bytes(meta_point) + sdsl::size_in_bytes(meta_rank) + 
                               sdsl::size_in_bytes(data_typ) + sdsl::size_in_bytes(data_val);
    }
};

}  // namespace h0