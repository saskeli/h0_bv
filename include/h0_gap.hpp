#include <immintrin.h>

#include <array>
#include <cstdint>
#include <memory>
#include <type_traits>

#include "internal.hpp"
#include "sdsl/int_vector.hpp"

namespace h0 {
namespace internal {

template <uint16_t block_size, uint16_t gap = block_size, bool flippable = true>
class mults {
 public:
  static_assert(gap > 0);
  static_assert(block_size < 32);
  typedef std::conditional<(block_size > 8), uint16_t, uint8_t>::type s_type;
  typedef std::conditional<(block_size > 16), uint32_t, s_type>::type dtype;

  const static constexpr auto bins = internal::binoms<block_size>();

  static constexpr uint64_t gap_elems() {
    uint64_t acc = 0;
    uint16_t types = flippable ? (block_size + 1) / 2 : block_size;
    for (uint16_t i = 0; i <= types; ++i) {
      acc += (bins[i] + gap - 1) / gap;
    }
    return acc;
  }

  const static constexpr std::array<
      uint64_t, 1 + (flippable ? (block_size + 1) / 2 : block_size)>
  goff() {
    std::array<uint64_t, 1 + (flippable ? (block_size + 1) / 2 : block_size)>
        ret;
    ret[0] = 0;
    for (uint16_t i = 1; i < ret.size(); ++i) {
      ret[i] = ret[i - 1] + (bins[i - 1] + gap - 1) / gap;
    }
    return ret;
  }

  const static constexpr std::array<uint64_t,
                                    (gap_elems() * block_size + 63) / 64>
  comp_repr_arr() {
    std::array<uint64_t, (gap_elems() * block_size + 63) / 64> ret;
    std::fill(ret.begin(), ret.end(), 0);
    ret[0] = 0;
    uint32_t idx = 1;
    for (uint16_t C = 1; C <= (flippable ? (block_size + 1) / 2 : block_size);
         ++C) {
      uint32_t v = (uint32_t(1) << C) - 1;
      uint32_t elems = (bins[C] + gap - 1) / gap;
      uint64_t i = idx * block_size;
      uint64_t off = i % 64;
      i /= 64;
      ret[i] |= uint64_t(v) << off;
      if (off + block_size > 64) {
        ret[i + 1] = uint64_t(v) >> (64 - off);
      }
      idx++;
      for (uint32_t ge = 1; ge < elems; ++ge) {
        for (uint32_t gi = 0; gi < gap; ++gi) {
          v = next(v);
        }
        i = idx * block_size;
        off = i % 64;
        i /= 64;
        ret[i] |= uint64_t(v) << off;
        if (off + block_size > 64) {
          ret[i + 1] = uint64_t(v) >> (64 - off);
        }
        idx++;
      }
    }
    return ret;
  }

  const static constexpr auto offsets = goff();
  const static constexpr auto gap_blocks = comp_repr_arr();

  static constexpr uint64_t get(uint64_t i) {
    const constexpr uint64_t mask = (uint64_t(1) << block_size) - 1;
    i *= block_size;
    uint64_t off = i % 64;
    i /= 64;
    uint64_t v = gap_blocks[i];
    v >>= off;
    if (block_size + off > 64) {
      uint64_t add = gap_blocks[i + 1];
      v |= add << (64 - off);
    }
    return v & mask;
  }

  template <bool blah = false>
  static constexpr uint64_t decode(uint16_t C, uint64_t offset) {
    const constexpr uint64_t mask = (uint64_t(1) << block_size) - 1;
    bool flip = false;
    if constexpr (flippable) {
      flip = C > (block_size + 1) / 2;
      C = flip ? block_size - C : C;
    }
    uint64_t idx = offset / gap;
    offset %= gap;
    uint64_t v = get(offsets[C] + idx);
    for (uint16_t i = 0; i < gap - 1; ++i) {
      if (i == offset) {
        break;
      }
      v = next(v);
    }
    if constexpr (flippable) {
      v = flip ? ~v : v;
      return v & mask;
    } else {
      return v;
    }
  }

  static std::unique_ptr<dtype[]> encode_array() {
    auto ret = std::make_unique<dtype[]>(uint64_t(1) << block_size);
    for (uint16_t C = 0; C <= block_size; ++C) {
      for (uint32_t off = 0; off < bins[C]; ++off) {
        ret[decode(C, off)] = off;
      }
    }
    return ret;
  }
};

template <uint16_t b>
const constexpr std::array<uint8_t, b + 1> f_widths() {
  std::array<uint8_t, b + 1> ret;
  auto bins = binoms<b>();
  for (uint16_t i = 0; i <= b; ++i) {
    ret[i] = 64 - __builtin_clzll(bins[i] - 1);
  }
  return ret;
}

}  // namespace internal

template <uint16_t block_size, uint16_t gap_size = block_size, uint16_t k = 32,
          bool flippable = true>
class h0_gap {
 private:
  static const constexpr internal::mults<block_size, gap_size> coder =
      internal::mults<block_size, gap_size>();
  static const constexpr auto widths = internal::f_widths<block_size>();
  static const constexpr uint16_t k_width = 64 - __builtin_clzll(block_size);
  sdsl::int_vector<> meta_point;
  sdsl::int_vector<> meta_rank;
  sdsl::int_vector<> data_typ;
  sdsl::bit_vector data_val;
  uint64_t size_;
  uint64_t sum_;

 public:
  h0_gap(const sdsl::bit_vector& bv) : sum_(0) {
    auto enc_arr = coder.encode_array();

    size_ = bv.size();
    uint64_t blocks = (size_ + block_size - 1) / block_size;
    uint64_t block = 0;
    uint64_t offset = 0;
    while (blocks--) {
      uint64_t elem = bv.get_int(block_size * block++, block_size);
      uint16_t p = __builtin_popcountll(elem);
      sum_ += p;
      offset += widths[p];
    }

    uint64_t v_siz = (size_ + k * block_size - 1) / (k * block_size);
    data_typ =
        sdsl::int_vector<>((size_ + block_size - 1) / block_size, 0, k_width);
    data_val = sdsl::bit_vector(offset);
    meta_point = sdsl::int_vector<>(v_siz, 0, 64 - __builtin_clzll(size_));
    meta_rank = sdsl::int_vector<>(v_siz, 0, 64 - __builtin_clzll(sum_));

    block = 0;
    sum_ = 0;
    offset = 0;
    blocks = (size_ + block_size - 1) / block_size;
    uint64_t s_block = 0;
    while (blocks--) {
      if (block % k == 0) {
        assert(meta_point.size() > s_block);
        meta_point[s_block] = offset;
        assert(meta_rank.size() > s_block);
        meta_rank[s_block++] = sum_;
      }
      uint64_t elem = bv.get_int(block_size * block, block_size);
      uint16_t pop = __builtin_popcountll(elem);
      uint64_t v = enc_arr[elem];
      // std::cerr << std::bitset<64>(elem) << " -> " << enc.first << ", "
      // << enc.second << std::endl;
      sum_ += pop;
      assert(data_typ.size() > block);
      data_typ[block++] = pop;
      data_val.set_int(offset, v, widths[pop]);
      // std::cerr << "                -> " << widths[enc.first] <<
      // std::endl;
      offset += widths[pop];
    }
  }

  bool access(uint64_t i) const {
    uint64_t block = i / block_size;
    auto bk = data_typ[block];
    if (bk == 0 || bk == block_size) {
      return bk != 0;
    }
    uint64_t s_block = block / k;
    uint16_t b_offset = block % k;
    i %= block_size;
    uint64_t offset = meta_point[s_block];
    for (uint64_t j = 0; j < k; ++j) {
      bk = data_typ[s_block * k + j];
      if (b_offset == j) {
        break;
      }
      offset += widths[bk];
    }
    auto f = data_val.get_int(offset, widths[bk]);
    return (coder.decode(bk, f) >> i) & 1;
  }

  uint64_t rank(uint64_t i) const {
    if (i >= size_) [[unlikely]] {
      return sum_;
    }
    uint64_t block = i / block_size;
    uint64_t s_block = block / k;
    uint16_t b_offset = block % k;
    i %= block_size;
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
    if (bk == 0 || bk == block_size) {
      return prefix_rank + i * (bk != 0);
    }
    auto f = data_val.get_int(offset, widths[bk]);
    uint64_t v = coder.decode(bk, f);
    v &= (uint64_t(1) << i) - 1;
    return prefix_rank + __builtin_popcountll(v);
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
    uint64_t res = k * a * block_size;
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
      res += block_size;
      s -= bk;
    }
    uint64_t v = coder.decode(bk, f);
    uint64_t pos = uint64_t(1) << (s - 1);
    return res + __builtin_ctzll(_pdep_u64(pos, v));
  }

  uint64_t size() const { return size_; }

  uint64_t bytes_size() const {
    return sizeof(h0_gap) + sdsl::size_in_bytes(meta_point) +
           sdsl::size_in_bytes(meta_rank) + sdsl::size_in_bytes(data_typ) +
           sdsl::size_in_bytes(data_val);
  }

  uint64_t C_bytes() const { return sdsl::size_in_bytes(data_typ); }

  uint64_t F_bytes() const { return sdsl::size_in_bytes(data_val); }

  uint64_t partial_sums() const {
    return sdsl::size_in_bytes(meta_point) + sdsl::size_in_bytes(meta_rank);
  }

  uint64_t lookup_tables() const {
    uint64_t ret = sizeof(decltype(coder)::bins) + sizeof(decltype(coder)::offsets) + sizeof(decltype(coder)::gap_blocks);
    return ret;
  }

  uint64_t object_bytes() const { return sizeof(*this); }
};

}  // namespace h0
