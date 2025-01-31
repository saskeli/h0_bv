#include <immintrin.h>

#include <algorithm>
#include <array>
#include <bitset>
#include <cstdint>
#include <iostream>
#include <vector>

#include "internal.hpp"
#include "sdsl/int_vector.hpp"

namespace h0 {
namespace internal {
class mults {
 public:
  static const constexpr uint16_t n = 64;
  static const constexpr std::array<uint64_t, n + 1> b64 = binoms<64>();
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
      tot += fp * b8[ks];
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
#ifndef RRR_NO_OPT
    if (f == 0) {
      if constexpr (b == 64) {
        if (k == 64) {
          return ~uint64_t(0);
        }
      }
      return (uint64_t(1) << k) - 1;
    }
    if constexpr (b == 64) {
      if (f == b64[k] - 1) {
        return ~((uint64_t(1) << (64 - k)) - 1);
      }
    }
    if constexpr (b == 32) {
      if (f == b32[k] - 1) {
        return ~((uint32_t(1) << (32 - k)) - 1);
      }
    }
    if constexpr (b == 16) {
      if (f == b16[k] - 1) {
        uint16_t ret = ~((uint16_t(1) << (16 - k)) - 1);
        // std::cerr << std::bitset<64>(ret) << " b16 last k = " << k <<
        // std::endl;
        return ret;
      }
    }
    if (k == 1) {
      return uint64_t(1) << f;
    }
    if (k == b - 1) {
      if (b == 64) {
        return (~uint64_t(0)) ^ (uint64_t(1) << (b - f - 1));
      } else {
        return ((uint64_t(1) << b) - 1) ^ (uint64_t(1) << (b - f - 1));
      }
    }
#endif
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
#ifndef RRR_NO_OPT
    if (i >= b) [[unlikely]] {
      return k;
    } else if (k == 0) [[unlikely]] {
      return 0;
    } else if (k == b) [[unlikely]] {
      return i;
    } else if (k == 1) [[unlikely]] {
      return f < i;
    } else if (k == b - 1) [[unlikely]] {
      return i - uint64_t((b - f - 1) < i);
    }
#endif
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
    if (i == b / 2) [[unlikely]] {
      return ks;
    }
    bool r = i > b / 2;
    bool ri = !r;
    k = r * kp + ri * ks;
    f = r * f / md + ri * f % md;
    i -= r * b / 2;
    return ks * r + rank<b / 2>(k, f, i);
  }

  template <uint16_t b>
  inline static constexpr uint64_t select(uint16_t k, uint64_t f, uint16_t s) {
    if constexpr (b == 8) {
      uint32_t by = byte_mapping<>::f_byte(k, f);
      uint16_t pos = uint16_t(1) << (s - 1);
      return __builtin_ctz(_pdep_u32(pos, by));
    }
#ifdef RRR_NO_OPT
    if (k == 0) [[unlikely]] {
      return 0;
    } else if (k == b) [[unlikely]] {
      return s - 1;
    } else if (k == 1) [[unlikely]] {
      return f;
    } else if (k == b - 1) [[unlikely]] {
      return s - 1 + ((b - f - 1) < s);
    }
#endif
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
    bool r = ks < s;
    bool ri = !r;
    k = r * kp + ri * ks;
    f = r * f / md + ri * f % md;
    return r * b / 2 + select<b / 2>(k, f, s - r * ks);
  }

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
    // std::cerr << "access(" << k << ", " << f << ", " << i << ")" <<
    // std::endl;
    if (f == 0) [[unlikely]] {
      return k > i;
    }
    if (f == b64[k] - 1) [[unlikely]] {
      /*std::cerr << "64: i = " << i << ", k = " << k << ", ";
      std::cerr << "f = " << b64[k] << " -> " << (i >= (n - k)) << std::endl;*/
      return i >= (n - k);
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
    if (f == 0) [[unlikely]] {
      return k > i;
    }
    if (f == b32[k] - 1) [[unlikely]] {
      return i >= (n / 2 - k);
    }
    if (k == 1) [[unlikely]] {
      return f == i;
    }
    if (k == n / 2 - 1) [[unlikely]] {
      return (31 - f) != i;
    }
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
    if (f == 0) [[unlikely]] {
      return k > i;
    }
    if (f == b16[k] - 1) [[unlikely]] {
      return i >= (n / 4 - k);
    }
    if (k == 1) [[unlikely]] {
      return f == i;
    }
    if (k == n / 4 - 1) [[unlikely]] {
      return (15 - f) != i;
    }
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

  inline static constexpr uint64_t select(uint16_t k, uint64_t f, uint16_t s) {
    return select<64>(k, f, s);
  }
};

template <uint16_t b>
const constexpr std::array<uint8_t, b + 1> f_widhts() {
  std::array<uint8_t, b + 1> ret;
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
    data_typ = sdsl::int_vector<>(
        size_ / block_width + (size_ % block_width > 0), 0, k_width);
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
    auto bk = data_typ[block];
    if (bk == 0 || bk == block_width) {
      return bk != 0;
    }
    uint64_t s_block = block / k;
    uint16_t b_offset = block % k;
    i %= block_width;
    uint64_t offset = meta_point[s_block];
    for (uint64_t j = 0; j < k; ++j) {
      bk = data_typ[s_block * k + j];
      if (b_offset == j) {
        break;
      }
      offset += widths[bk];
    }
    auto f = data_val.get_int(offset, widths[bk]);
    return coder.access(bk, f, i);
  }

  uint64_t rank(uint64_t i) const {
    if (i >= size_) [[unlikely]] {
      return sum_;
    }
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
    return sizeof(h0_bv) + sdsl::size_in_bytes(meta_point) +
           sdsl::size_in_bytes(meta_rank) + sdsl::size_in_bytes(data_typ) +
           sdsl::size_in_bytes(data_val);
  }

  uint64_t C_bytes() const { return sdsl::size_in_bytes(data_typ); }

  uint64_t F_bytes() const { return sdsl::size_in_bytes(data_val); }

  uint64_t partial_sums() const {
    return sdsl::size_in_bytes(meta_point) + sdsl::size_in_bytes(meta_rank);
  }

  uint64_t lookup_tables() const {
    uint64_t ret = sizeof(internal::mults::b64) + sizeof(internal::mults::b32) +
                   sizeof(internal::mults::b16) + sizeof(internal::mults::b8) +
                   sizeof(internal::mults::f_lim64) +
                   sizeof(internal::mults::f_lim32) +
                   sizeof(internal::mults::f_lim16) + sizeof(widths);
    return ret;
  }

  uint64_t object_bytes() const { return sizeof(*this); }
};

}  // namespace h0