from math import comb, ceil
import sys

def is_necklace(s):
    if s[-1] == '0':
        return False
    if s[0] == '1':
        return False
    lrl = 0
    while s[lrl] == '0':
        lrl += 1
    larl = 0
    lars = []
    idx = lrl
    while idx < len(s):
        nrl = 0
        sidx = idx
        while s[idx] == '0':
            nrl += 1
            idx += 1
        if nrl > larl:
            larl = nrl
            lars = [sidx]
        elif nrl == lrl:
            lars.append(sidx)
        idx += 1
    if lrl < larl:
        return False
    if lrl > larl:
        return True
    for i in range(len(s)):
        nlars = []
        if s[(lrl + i) % len(s)] == '0':
            for rs in lars:
                if s[(rs + lrl + i) % len(s)] == '0':
                    nlars.append(rs)
        else:
            for rs in lars:
                if s[(rs + lrl + i) % len(s)] == '0':
                    return False
                nlars.append(rs)
        if len(nlars) == 0:
            return True
        lars = nlars
    return True

def shift(s, a, b):
    sp = s[:b - 1] + s[a - 1] + s[b - 1:a - 1] + s[a:]
    return sp

def next_s(s):
    idx = 0
    while s[idx] == '0':
        idx += 1
    zeros = idx
    while idx < len(s) and s[idx] == '1':
        idx += 1
    ones = idx - zeros
    gamma = s[ones + zeros:]
    if len(gamma) == 0:
        return shift(s, zeros + ones, 1 + (zeros + 1) // 2)
    beta = shift(s, zeros + ones + 2, zeros + ones + 1)
    if s[ones + zeros + 1] == '0' or is_necklace(beta) == False:
        return shift(s, zeros + ones + 1, 1)
    i = 1
    while i < zeros + ones + 2:
        thing = shift(s, zeros + ones + 2, i)
        if is_necklace(thing):
            break
        i += 1
    return shift(s, zeros + ones + 2, i)

def aperiodic_prefix(s):
    for i in range(1, len(s)):
        if s == s[i:] + s[:i]:
            return s[:i]
    return s

def weighted_deBrujin_seq(starts):
    s = starts
    neck = []
    while True:
        s = next_s(s)
        neck.append(s)
        if s == starts:
            break
    fddbs = ""
    for s in reversed(neck):
        fddbs += aperiodic_prefix(s)
    return fddbs

def cpp_header(n):
    print("""#include <array>
#include <cstdint>
#include <endian.h>
#include <immintrin.h>
#include <sdsl/int_vector.hpp>

#include "internal.hpp"

namespace h0 {
namespace internal {

class wdbs {""")
    max_pop = (n + 1) // 2
    arr_size = ceil((sum(comb(n, i) + 1 for i in range(1, max_pop + 1)) + n - 1) / 8)
    
    print(f"    static const constexpr std::array<uint8_t, {arr_size}> data = {'{'}", end="")
    arridx = 0
    offsets = [0, 0]

    one_seq = "0" * (n - 1) + "1"
    one_blocks = ceil(len(one_seq) / 8)
    partial_block = ""
    for j in range(one_blocks):
        if arridx % 12 == 0:
            print("\n        ", end="")
        block = one_seq[8 * j:8 * j + 8]
        if j == one_blocks - 1 and len(block) < 8:
            partial_block = block
            break
        print(f"{int(block, 2)}u, ", end="")
        arridx += 1
        
    for i in range(2, max_pop + 1):
        ins = weighted_deBrujin_seq("0" * (n - i) + "1" * i)
        ins = partial_block + ins + ("0" if i < max_pop else "0" * (n - max_pop) + "1" * max_pop)
        offsets.append(arridx * 8 + len(partial_block))
        partial_block = ""
        blocks = ceil(len(ins) / 8)
        for j in range(blocks):
            if arridx % 12 == 0:
                print("\n        ", end="")
            block = ins[8 * j:8 * j + 8]
            if j == blocks - 1 and len(block) < 8:
                partial_block = block
                break
            print(f"{int(block, 2)}u", end="")
            arridx += 1
            if j < blocks - 1 or i < max_pop:
                print(", ", end="")
    if len(partial_block) > 0:
        partial_block += "0" * (8 - len(partial_block))
        print(f"{int(partial_block, 2)}u", end="")
    print("};")
    
    i_typ = "uint8_t"
    if offsets[-1] >= 2**32:
        i_typ = "uint64_t"
    if offsets[-1] >= 2**16:
        i_typ = "uint32_t"
    elif offsets[-1] >= 2**8:
        i_typ = "uint16_t"

    d_typ = 8
    if n > 32:
        d_typ = 64
    elif n > 16:
        d_typ = 32
    elif n > 8:
        d_typ = 16

    print(f"    static const constexpr std::array<{i_typ}, {len(offsets)}> offs = {'{'}{', '.join(str(v) for v in offsets)}{'}'};")

    print("   public:")
    print(f"    static constexpr uint{d_typ}_t decode(uint16_t C, {i_typ} off) {'{'}")
    print(f"        const constexpr uint{d_typ}_t mask = (uint{d_typ}_t(1) << {n}) - 1;")
    print(f"        bool flip = C > {ceil(n / 2)};")
    print(f"        C = flip ? {n} - C : C;")
    print("        off += offs[C];")
    print("        uint32_t idx = off / 8;")
    print("        off %= 8;")
    print("        uint64_t w = be64toh(reinterpret_cast<const uint64_t*>(data.data() + idx)[0]);")
    print(f"        w = (w >> ({64 - n} - off)) & mask;")
    print("        w ^= (C - __builtin_popcountll(w)) & uint64_t(1);")
    print("        w = flip ? ~w & mask: w;")
    print("        return w;")
    print("    }\n")
    print("};")

    i_typ = "uint8_t"
    if n >= 32:
        i_typ = "uint64_t"
    elif n >= 16:
        i_typ = "uint32_t"
    elif n >= 8:
        i_typ = "uint16_t"
    print(f"std::unique_ptr<uint{d_typ}_t[]>  encode_array() {'{'}")
    print(f"    auto bin = binoms<{n}>();")
    print(f"    auto ret = std::make_unique<uint{d_typ}_t[]>({2**n});")
    print(f"    for (uint16_t C = 0; C <= {n}; ++C) {'{'}")
    print(f"        for ({i_typ} off = 0; off < bin[C]; ++off) {'{'}")
    print("            ret[wdbs::decode(C, off)] = off;")
    print("        }")
    print("    }")
    print("    return ret;\n}")
    print(f"""
const constexpr std::array<uint8_t, {n + 1}> f_widths() {'{'}
    std::array<uint8_t, {n + 1}> ret;
    auto bins = binoms<{n}>();
    for (uint16_t i = 0; i <= {n}; ++i) {'{'}
        ret[i] = 64 - __builtin_clzll(bins[i] - 1);
    {'}'}
    return ret;
{'}'}

{'}'} // namespace internal

template <uint16_t k = 32>
class h0_wdb {'{'}
   private:
    static const constexpr internal::wdbs coder = internal::wdbs();
    static const constexpr auto widths = internal::f_widths();
    static const constexpr uint16_t k_width = 64 - __builtin_clzll({n});
    sdsl::int_vector<> meta_point;
    sdsl::int_vector<> meta_rank;
    sdsl::int_vector<> data_typ;
    sdsl::bit_vector data_val;
    uint64_t size_;
    uint64_t sum_;

   public:
    h0_wdb(const sdsl::bit_vector& bv) : sum_(0) {'{'}
        auto enc_arr = internal::encode_array();

        size_ = bv.size();
        uint64_t blocks = (size_ + {n} - 1) / {n};
        uint64_t block = 0;
        uint64_t offset = 0;
        while (blocks--) {'{'}
            uint64_t elem = bv.get_int({n} * block++, {n});
            uint16_t p = __builtin_popcountll(elem);
            sum_ += p;
            offset += widths[p];
        {'}'}

        uint64_t v_siz = (size_ + k * {n} - 1) / (k * {n});
        data_typ = sdsl::int_vector<>((size_ + {n} - 1) / {n}, 0, k_width);
        data_val = sdsl::bit_vector(offset);
        meta_point = sdsl::int_vector<>(v_siz, 0, 64 - __builtin_clzll(size_));
        meta_rank = sdsl::int_vector<>(v_siz, 0, 64 - __builtin_clzll(sum_));

        block = 0;
        sum_ = 0;
        offset = 0;
        blocks = (size_ + {n} - 1) / {n};
        uint64_t s_block = 0;
        while (blocks--) {'{'}
            if (block % k == 0) {'{'}
                assert(meta_point.size() > s_block);
                meta_point[s_block] = offset;
                assert(meta_rank.size() > s_block);
                meta_rank[s_block++] = sum_;
            {'}'}
            uint64_t elem = bv.get_int({n} * block, {n});
            uint16_t pop = __builtin_popcountll(elem);
            uint64_t v = enc_arr[elem];
            sum_ += pop;
            assert(data_typ.size() > block);
            data_typ[block++] = pop;
            data_val.set_int(offset, v, widths[pop]);
            offset += widths[pop];
        {'}'}
    {'}'}

    bool access(uint64_t i) const {'{'}
        uint64_t block = i / {n};
        uint64_t s_block = block / k;
        uint16_t b_offset = block % k;
        i %= {n};
        uint64_t offset = meta_point[s_block];
        for (uint64_t j = 0; j < k; ++j) {'{'}
            if (b_offset == j) {'{'}
                break;
            {'}'}
            auto bk = data_typ[s_block * k + j];
            offset += widths[bk];
        {'}'}
        auto bk = data_typ[block];
        auto f = data_val.get_int(offset, widths[bk]);
        return (coder.decode(bk, f) >> i) & 1;
    {'}'}

    uint64_t rank(uint64_t i) const {'{'}
        if (i >= size_) [[unlikely]] {'{'}
            return sum_;
        {'}'}
        uint64_t block = i / {n};
        uint64_t s_block = block / k;
        uint16_t b_offset = block % k;
        i %= {n};
        uint64_t offset = meta_point[s_block];
        uint64_t prefix_rank = meta_rank[s_block];
        for (uint64_t j = 0; j < k; ++j) {'{'}
            if (b_offset == j) {'{'}
                break;
            {'}'}
            auto bk = data_typ[s_block * k + j];
            prefix_rank += bk;
            offset += widths[bk];
        {'}'}
        auto bk = data_typ[block];
        auto f = data_val.get_int(offset, widths[bk]);
        uint64_t v = coder.decode(bk, f);
        v &= (uint64_t(1) << i) - 1;
        return prefix_rank + __builtin_popcountll(v);
    {'}'}

    uint64_t select(uint64_t s) const {'{'}
        uint64_t a = 0;
        uint64_t b = meta_point.size() - 1;
        while (a < b) {'{'}
            uint64_t m = (a + b + 1) / 2;
            if (meta_rank[m] >= s) {'{'}
                b = m - 1;
            {'}'} else {'{'}
                a = m;
            {'}'}
        {'}'}
        uint64_t offset = meta_point[a];
        s -= meta_rank[a];
        uint64_t res = k * a * {n};
        uint16_t bk = 0;
        uint64_t f = 0;
        uint64_t block = a * k;
        for (uint16_t i = 0; i < k; ++i) {'{'}
            bk = data_typ[block++];
            if (bk >= s) {'{'}
                f = data_val.get_int(offset, widths[bk]);
                break;
            {'}'}
            offset += widths[bk];
            res += {n};
            s -= bk;
        {'}'}
        uint64_t v = coder.decode(bk, f);
        uint64_t pos = uint64_t(1) << (s - 1);
        return res + __builtin_ctzll(_pdep_u64(pos, v));
    {'}'}

    uint64_t size() const {'{'} return size_; {'}'}

    uint64_t bytes_size() const {'{'}
        return sizeof(h0_wdb) + sdsl::size_in_bytes(meta_point) + sdsl::size_in_bytes(meta_rank) + 
                               sdsl::size_in_bytes(data_typ) + sdsl::size_in_bytes(data_val);
    {'}'}
{'}'};

{'}'}  // namespace h0
""")

def full_header(n):
    print("""#include <array>
#include <cstdint>
#include <endian.h>
#include <immintrin.h>
#include <sdsl/int_vector.hpp>

#include "internal.hpp"

namespace h0 {
namespace internal {

class wdbs {""")
    arr_size = ceil((sum(comb(n, i) + 1 for i in range(1, n)) + n) / 8)
    
    print(f"    static const constexpr std::array<uint8_t, {arr_size}> data = {'{'}", end="")
    arridx = 0
    offsets = [0, 0]

    one_seq = "0" * (n - 1) + "1"
    one_blocks = ceil(len(one_seq) / 8)
    partial_block = ""
    for j in range(one_blocks):
        if arridx % 12 == 0:
            print("\n        ", end="")
        block = one_seq[8 * j:8 * j + 8]
        if j == one_blocks - 1 and len(block) < 8:
            partial_block = block
            break
        print(f"{int(block, 2)}u, ", end="")
        arridx += 1
        
    for i in range(2, n):
        ins = weighted_deBrujin_seq("0" * (n - i) + "1" * i)
        ins = partial_block + ins + "0"
        offsets.append(arridx * 8 + len(partial_block))
        partial_block = ""
        if i == n - 1:
            offsets.append(arridx * 8 + len(ins))
            ins += "1" * n
        blocks = ceil(len(ins) / 8)
        for j in range(blocks):
            if arridx % 12 == 0:
                print("\n        ", end="")
            block = ins[8 * j:8 * j + 8]
            if j == blocks - 1 and len(block) < 8:
                partial_block = block
                break
            print(f"{int(block, 2)}u", end="")
            arridx += 1
            if j < blocks - 1 or i < n:
                print(", ", end="")
    if len(partial_block) > 0:
        partial_block += "0" * (8 - len(partial_block))
        print(f"{int(partial_block, 2)}u", end="")
    print("};")
    
    i_typ = "uint8_t"
    if offsets[-1] >= 2**32:
        i_typ = "uint64_t"
    if offsets[-1] >= 2**16:
        i_typ = "uint32_t"
    elif offsets[-1] >= 2**8:
        i_typ = "uint16_t"

    d_typ = 8
    if n > 32:
        d_typ = 64
    elif n > 16:
        d_typ = 32
    elif n > 8:
        d_typ = 16

    print(f"    static const constexpr std::array<{i_typ}, {len(offsets)}> offs = {'{'}{', '.join(str(v) for v in offsets)}{'}'};")

    print("   public:")
    print(f"    static constexpr uint{d_typ}_t decode(uint16_t C, {i_typ} off) {'{'}")
    print(f"        const constexpr uint{d_typ}_t mask = (uint{d_typ}_t(1) << {n}) - 1;")
    print("        off += offs[C];")
    print("        uint32_t idx = off / 8;")
    print("        off %= 8;")
    print("        uint64_t w = be64toh(reinterpret_cast<const uint64_t*>(data.data() + idx)[0]);")
    print(f"        w = (w >> ({64 - n} - off)) & mask;")
    print("        w ^= (C - __builtin_popcountll(w)) & uint64_t(1);")
    print("        return w;")
    print("    }\n")
    print("};")

    i_typ = "uint8_t"
    if n >= 32:
        i_typ = "uint64_t"
    elif n >= 16:
        i_typ = "uint32_t"
    elif n >= 8:
        i_typ = "uint16_t"
    print(f"std::unique_ptr<uint{d_typ}_t[]>  encode_array() {'{'}")
    print(f"    auto bin = binoms<{n}>();")
    print(f"    auto ret = std::make_unique<uint{d_typ}_t[]>({2**n});")
    print(f"    for (uint16_t C = 0; C <= {n}; ++C) {'{'}")
    print(f"        for ({i_typ} off = 0; off < bin[C]; ++off) {'{'}")
    print("            ret[wdbs::decode(C, off)] = off;")
    print("        }")
    print("    }")
    print("    return ret;\n}")
    print(f"""
const constexpr std::array<uint8_t, {n + 1}> f_widths() {'{'}
    std::array<uint8_t, {n + 1}> ret;
    auto bins = binoms<{n}>();
    for (uint16_t i = 0; i <= {n}; ++i) {'{'}
        ret[i] = 64 - __builtin_clzll(bins[i] - 1);
    {'}'}
    return ret;
{'}'}

{'}'} // namespace internal

template <uint16_t k = 32>
class h0_wdb {'{'}
   private:
    static const constexpr internal::wdbs coder = internal::wdbs();
    static const constexpr auto widths = internal::f_widths();
    static const constexpr uint16_t k_width = 64 - __builtin_clzll({n});
    sdsl::int_vector<> meta_point;
    sdsl::int_vector<> meta_rank;
    sdsl::int_vector<> data_typ;
    sdsl::bit_vector data_val;
    uint64_t size_;
    uint64_t sum_;

   public:
    h0_wdb(const sdsl::bit_vector& bv) : sum_(0) {'{'}
        auto enc_arr = internal::encode_array();

        size_ = bv.size();
        uint64_t blocks = (size_ + {n} - 1) / {n};
        uint64_t block = 0;
        uint64_t offset = 0;
        while (blocks--) {'{'}
            uint64_t elem = bv.get_int({n} * block++, {n});
            uint16_t p = __builtin_popcountll(elem);
            sum_ += p;
            offset += widths[p];
        {'}'}

        uint64_t v_siz = (size_ + k * {n} - 1) / (k * {n});
        data_typ = sdsl::int_vector<>((size_ + {n} - 1) / {n}, 0, k_width);
        data_val = sdsl::bit_vector(offset);
        meta_point = sdsl::int_vector<>(v_siz, 0, 64 - __builtin_clzll(size_));
        meta_rank = sdsl::int_vector<>(v_siz, 0, 64 - __builtin_clzll(sum_));

        block = 0;
        sum_ = 0;
        offset = 0;
        blocks = (size_ + {n} - 1) / {n};
        uint64_t s_block = 0;
        while (blocks--) {'{'}
            if (block % k == 0) {'{'}
                assert(meta_point.size() > s_block);
                meta_point[s_block] = offset;
                assert(meta_rank.size() > s_block);
                meta_rank[s_block++] = sum_;
            {'}'}
            uint64_t elem = bv.get_int({n} * block, {n});
            uint16_t pop = __builtin_popcountll(elem);
            uint64_t v = enc_arr[elem];
            sum_ += pop;
            assert(data_typ.size() > block);
            data_typ[block++] = pop;
            data_val.set_int(offset, v, widths[pop]);
            offset += widths[pop];
        {'}'}
    {'}'}

    bool access(uint64_t i) const {'{'}
        uint64_t block = i / {n};
        uint64_t s_block = block / k;
        uint16_t b_offset = block % k;
        i %= {n};
        uint64_t offset = meta_point[s_block];
        for (uint64_t j = 0; j < k; ++j) {'{'}
            if (b_offset == j) {'{'}
                break;
            {'}'}
            auto bk = data_typ[s_block * k + j];
            offset += widths[bk];
        {'}'}
        auto bk = data_typ[block];
        auto f = data_val.get_int(offset, widths[bk]);
        return (coder.decode(bk, f) >> i) & 1;
    {'}'}

    uint64_t rank(uint64_t i) const {'{'}
        if (i >= size_) [[unlikely]] {'{'}
            return sum_;
        {'}'}
        uint64_t block = i / {n};
        uint64_t s_block = block / k;
        uint16_t b_offset = block % k;
        i %= {n};
        uint64_t offset = meta_point[s_block];
        uint64_t prefix_rank = meta_rank[s_block];
        for (uint64_t j = 0; j < k; ++j) {'{'}
            if (b_offset == j) {'{'}
                break;
            {'}'}
            auto bk = data_typ[s_block * k + j];
            prefix_rank += bk;
            offset += widths[bk];
        {'}'}
        auto bk = data_typ[block];
        auto f = data_val.get_int(offset, widths[bk]);
        uint64_t v = coder.decode(bk, f);
        v &= (uint64_t(1) << i) - 1;
        return prefix_rank + __builtin_popcountll(v);
    {'}'}

    uint64_t select(uint64_t s) const {'{'}
        uint64_t a = 0;
        uint64_t b = meta_point.size() - 1;
        while (a < b) {'{'}
            uint64_t m = (a + b + 1) / 2;
            if (meta_rank[m] >= s) {'{'}
                b = m - 1;
            {'}'} else {'{'}
                a = m;
            {'}'}
        {'}'}
        uint64_t offset = meta_point[a];
        s -= meta_rank[a];
        uint64_t res = k * a * {n};
        uint16_t bk = 0;
        uint64_t f = 0;
        uint64_t block = a * k;
        for (uint16_t i = 0; i < k; ++i) {'{'}
            bk = data_typ[block++];
            if (bk >= s) {'{'}
                f = data_val.get_int(offset, widths[bk]);
                break;
            {'}'}
            offset += widths[bk];
            res += {n};
            s -= bk;
        {'}'}
        uint64_t v = coder.decode(bk, f);
        uint64_t pos = uint64_t(1) << (s - 1);
        return res + __builtin_ctzll(_pdep_u64(pos, v));
    {'}'}

    uint64_t size() const {'{'} return size_; {'}'}

    uint64_t bytes_size() const {'{'}
        return sizeof(h0_wdb) + sdsl::size_in_bytes(meta_point) + sdsl::size_in_bytes(meta_rank) + 
                               sdsl::size_in_bytes(data_typ) + sdsl::size_in_bytes(data_val);
    {'}'}
{'}'};

{'}'}  // namespace h0
""")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("n is required")
    else:
        if len(sys.argv) > 2:
            full_header(int(sys.argv[1]))
        else:
            cpp_header(int(sys.argv[1]))