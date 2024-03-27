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

#include "internal.hpp"

class wdbs {""")
    max_pop = n // 2
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
    print(f"    static const constexpr uint{d_typ}_t decode(uint16_t C, {i_typ} off) {'{'}")
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
    print(f"constexpr std::array<uint{d_typ}_t, {2**n}> encode_array() {'{'}")
    print(f"    auto bin = binoms<{n}>();")
    print(f"    std::array<uint{d_typ}_t, {2**n}> ret;")
    print(f"    for (uint16_t C = 0; C <= {n // 2}; ++C) {'{'}")
    print(f"        for ({i_typ} off = 0; off < bin[C]; ++off) {'{'}")
    print("            ret[wdbs::decode(C, off)] = off;")
    print("        }")
    print("    }")
    print(f"    for ({i_typ} i = 0; i < {2**n}; ++i) {'{'}")
    print(f"        if (__builtin_popcountll(i) > {n // 2}) {'{'}")
    print(f"            ret[i] = ret[(~i) & 0b{'1' * n}];")
    print("        }")
    print("    }")
    print("    return ret;\n}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("n is required")
    else:
        cpp_header(int(sys.argv[1]))