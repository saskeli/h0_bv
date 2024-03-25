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
    fddbs += aperiodic_prefix(neck[-1][:-1])
    return fddbs

def cpp_header(n):
    print("""#include <array>
#include <cstdint>

static const constexpr class {""")
    max_pop = n // 2
    arr_size = 0
    offsets = []
    for i in range(2, max_pop + 1):
        offsets.append(arr_size)
        arr_size += ceil((comb(n, i) + n - 1) / 8)
    
    print(f"    std::array<uint8_t, {arr_size}> data = {'{'}", end="")
    arridx = 0
    for i in range(2, max_pop + 1):
        ins = weighted_deBrujin_seq("0" * (n - i) + "1" * i)
        blocks = ceil(len(ins) / 8)
        for j in range(blocks):
            if arridx % 10 == 0:
                print("\n        ", end="")
            block = ins[8 * j:8 * j + 8]
            if j == blocks - 1:
                block += "0" * (8 - len(block))
            print(f"{int(block[::-1], 2)}u", end="")
            print("};\n" if i == max_pop and j == blocks - 1 else ", ", end="")
            arridx += 1
    
    print(f"    std::array<uint32_t, {len(offsets)}> offs = {'{'}{', '.join(str(v) for v in offsets)}{'}'};")

    print("   public:")
    print("    uint32_t decode(uint16_t C, uint32_t off) const {")
    print(f"        const constexpr uint32_t mask = (uint32_t(1) << {n}) - 1;")
    print("        uint32_t idx = off / 8 + offs[C - 2];")
    print("        off %= 8;")
    print("        uint64_t w = reinterpret_cast<const uint64_t*>(data.data() + idx)[0];")
    print("        w = (w >> off) & mask;")
    print("        w ^= (C - __builtin_popcountll(w)) & uint64_t(1);")
    print("        return w;")
    print("    }\n")
    print("    std::pair<uint8_t*, uint32_t> get_arr(uint16_t C) {")
    print("        C -= 2;")
    print("        uint8_t* ptr = data.data() + offs[C];")
    print("        uint32_t size = C == offs.size() - 2 ? data.size() : offs[C];")
    print("        size -= offs[C];")
    print("        return {ptr, size};")
    print("    }")
    print("} wbds;")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("n is required")
    else:
        cpp_header(int(sys.argv[1]))