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
    static const constexpr uint16_t n = 63;
    static const constexpr std::array<uint64_t, 64> b63 = binoms<63>();
    static const constexpr std::array<uint64_t, 57> b56 = binoms<56>();
    static const constexpr std::array<uint64_t, 49> b48 = binoms<48>();
    static const constexpr std::array<uint64_t, 41> b40 = binoms<40>();
    static const constexpr std::array<uint64_t, 33> b32 = binoms<32>();
    static const constexpr std::array<uint64_t, 25> b24 = binoms<24>();
    static const constexpr std::array<uint64_t, 17> b16 = binoms<16>();
    static const constexpr std::array<uint64_t, 9> b8 = binoms<8>();
    static const constexpr std::array<uint64_t, 8> b7 = binoms<7>();
    static const constexpr std::array<std::array<uint64_t, 9>, 64>
        f_lim63 = f_lims_it<63>();
    static const constexpr std::array<std::array<uint64_t, 9>, 57>
        f_lim56 = f_lims_it<56>();
    static const constexpr std::array<std::array<uint64_t, 9>, 49>
        f_lim48 = f_lims_it<48>();
    static const constexpr std::array<std::array<uint64_t, 9>, 41>
        f_lim40 = f_lims_it<40>();
    static const constexpr std::array<std::array<uint64_t, 9>, 33>
        f_lim32 = f_lims_it<32>();
    static const constexpr std::array<std::array<uint64_t, 9>, 25>
        f_lim24 = f_lims_it<24>();
    static const constexpr std::array<std::array<uint64_t, 9>, 17>
        f_lim16 = f_lims_it<16>();

    template <uint16_t b>
    inline static constexpr uint64_t encode(uint16_t k, uint64_t bv) {
        if constexpr (b <= 8) {
            //std::cerr << "byte f: " << std::bitset<8>(bv) << " -> " << byte_mapping<>::get_f(bv) << std::endl;
            return byte_mapping<>::get_f(bv);
        } else {
            if (k == 0) {
                //std::cerr << b << " 0f: " << 0 << std::endl;
                return 0;
            } else if (k == 1) {
                //std::cerr << b << " 1f: " << __builtin_ctzll(bv) << std::endl;
                return __builtin_ctzll(bv);
            } else if (k == b) {
                return 0;
            } else if (k == b - 1) {
                return b - __builtin_ctzll(~bv) - 1;
            }
            const constexpr uint16_t next = b == 63 ? 56 : b - 8;
            const constexpr uint64_t mask = (uint64_t(1) << next) - 1;
            uint64_t bvs = bv & mask;
            uint16_t ks = __builtin_popcountll(bvs);
            uint16_t kp = k - ks;
            //std::cerr << b << " bit f: " << std::bitset<8>(bv >> next) << " -> " << byte_mapping<>::get_f(bv >> next) << std::endl;
            uint64_t fp = byte_mapping<>::get_f(bv >> next);
            uint64_t fs = encode<next>(ks, bvs);
            //std::cerr << b << " suffi f: " << fs << std::endl;
            uint64_t tot = 0;
            if constexpr (b == 63) {
                tot += f_lim63[k][kp];
            } else if constexpr (b == 56) {
                tot += f_lim56[k][kp];
            } else if constexpr (b == 48) {
                tot += f_lim48[k][kp];
            } else if constexpr (b == 40) {
                tot += f_lim40[k][kp];
            } else if constexpr (b == 32) {
                tot += f_lim32[k][kp];
            } else if constexpr (b == 24) {
                tot += f_lim24[k][kp];
            } else {
                tot += f_lim16[k][kp];
            }
            /*if constexpr (b == 63) {
                std::cerr << "low k offset: " << tot << std::endl;
                std::cerr << "prefix contrib: " << fp * b56[ks] << std::endl;
                std::cerr << "suffix contrib: " << fs << std::endl;
            }//*/

            if constexpr (b == 63) {
                tot += fp * b56[ks];
            } else if constexpr (b == 56) {
                tot += fp * b48[ks];
            } else if constexpr (b == 48) {
                tot += fp * b40[ks];
            } else if constexpr (b == 40) {
                tot += fp * b32[ks];
            } else if constexpr (b == 32) {
                tot += fp * b24[ks];
            } else if constexpr (b == 24) {
                tot += fp * b16[ks];
            } else {
                tot += fp * b8[ks];
            }
            tot += fs;
            //std::cerr << b << " f: " << tot << std::endl;
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
        /*std::cerr << "f = " << f << " compares to ";
        for (auto v : fl) {
            std::cerr << v << ", ";
        } 
        std::cerr << "leading to kp = " << kp << std::endl;*/
        f -= f_lim;
    }

   public:
    inline static constexpr std::pair<uint16_t, uint64_t> encode(uint64_t bv) {
        /*std::cerr << "encoding " << bv << "(" << std::bitset<8>(bv >> 56) << " "
                  << std::bitset<8>(bv >> 48) << " " << std::bitset<8>(bv >> 40)
                  << " " << std::bitset<8>(bv >> 32) << " "
                  << std::bitset<8>(bv >> 24) << " " << std::bitset<8>(bv >> 16)
                  << " " << std::bitset<8>(bv >> 8) << " " << std::bitset<8>(bv)
                  << ")" << std::endl;//*/

        uint16_t k = __builtin_popcountll(bv);
        auto f = encode<63>(k, bv);
        //std::cerr << " -> " << k << ", " << f << std::endl;
        return {k, f};
    }
    

    inline static constexpr uint64_t decode(uint16_t k, uint64_t f) {
        //std::cout << "decode " << k << ", " << f << std::endl;
        uint16_t kp;
        uint64_t f_lim;
        const constexpr uint64_t mask63 = ~(uint64_t(1) << 63);

        if (f == 0) {
            return (uint64_t(1) << k) - 1;
        }
        if (f == b63[k] - 1) {
            uint64_t val = uint64_t(1) << (n - k);
            val -= 1;
            return (~val) & mask63;
        }
        if (k == 1) {
            return uint64_t(1) << f;
        }
        if (k == (n - 1)) {
            uint64_t val = uint64_t(1) << (62 - f);
            val = ~val;
            return val & mask63;
        }

        bs(f_lim63[k], kp, f_lim, f);
        k -= kp;
        uint64_t md = b56[k];
        //std::cout << "63 flim: " << f_lim << ", pf: " << f / md << ", sf: " << f % md << ", kp: " << kp << std::endl;
        uint64_t res = byte_mapping<>::f_byte(kp, f / md);
        //std::cout << " res: " << std::bitset<64>(res) << std::endl;
        f %= md;

        bs(f_lim56[k], kp, f_lim, f);
        k -= kp;
        md = b48[k];
        //std::cout << "56 pf: " << f / md << ", sf: " << f % md << std::endl;
        res = (res << 8) | byte_mapping<>::f_byte(kp, f / md);
        //std::cout << " res: " << std::bitset<64>(res) << std::endl;
        f %= md;

        bs(f_lim48[k], kp, f_lim, f);
        k -= kp;
        md = b40[k];
        //std::cout << "48 pf: " << f / md << ", sf: " << f % md << std::endl;
        res = (res << 8) | byte_mapping<>::f_byte(kp, f / md);
        //std::cout << " res: " << std::bitset<64>(res) << std::endl;
        f %= md;

        bs(f_lim40[k], kp, f_lim, f);
        k -= kp;
        md = b32[k];
        //std::cout << "40 pf: " << f / md << ", sf: " << f % md << std::endl;
        res = (res << 8) | byte_mapping<>::f_byte(kp, f / md);
        //std::cout << " res: " << std::bitset<64>(res) << std::endl;
        f %= md;

        bs(f_lim32[k], kp, f_lim, f);
        k -= kp;
        md = b24[k];
        //std::cout << "32 pf: " << f / md << ", sf: " << f % md << std::endl;
        res = (res << 8) | byte_mapping<>::f_byte(kp, f / md);
        //std::cout << " res: " << std::bitset<64>(res) << std::endl;
        f %= md;

        bs(f_lim24[k], kp, f_lim, f);
        k -= kp;
        md = b16[k];
        //std::cout << "24 pf: " << f / md << ", sf: " << f % md << std::endl;
        res = (res << 8) | byte_mapping<>::f_byte(kp, f / md);
        //std::cout << " res: " << std::bitset<64>(res) << std::endl;
        f %= md;

        bs(f_lim16[k], kp, f_lim, f);
        k -= kp;
        md = b8[k];
        //std::cout << "16 pf: " << f / md << ", sf: " << f % md << std::endl;
        res = (res << 8) | byte_mapping<>::f_byte(kp, f / md);
        //std::cout << " res: " << std::bitset<64>(res) << std::endl;

        res = (res << 8) | byte_mapping<>::f_byte(k, f % md);
        //std::cout << " res: " << std::bitset<64>(res) << std::endl;
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
    static const constexpr uint16_t block_width = 63;
    static const constexpr internal::mults coder = internal::mults();
    static const constexpr auto widths = internal::f_widhts<block_width>();
    static const constexpr uint16_t k_width = 6;
    sdsl::int_vector<> meta_point;
    sdsl::int_vector<> meta_rank;
    sdsl::int_vector<> data_typ;
    sdsl::bit_vector data_val;
    uint64_t size_;
    uint64_t sum_;

   public:
    h0_bv(const sdsl::bit_vector& bv) : sum_(0) {
        size_ = bv.size();
        uint64_t blocks = size_ / block_width + (size_ % block_width ? 1 : 0);
        uint64_t block = 0;
        uint64_t offset = 0;
        while (blocks--) {
            uint64_t elem = bv.get_int(block_width * block++, block_width);
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
        std::vector<uint64_t> metap_cont;
        std::vector<uint64_t> metar_cont;
        std::vector<uint16_t> datat_cont;
        std::vector<uint64_t> datav_cont;

        block = 0;
        sum_ = 0;
        offset = 0;
        blocks = size_ / block_width + (size_ % block_width ? 1 : 0);
        uint64_t s_block = 0;
        while (blocks--) {
            if (block % k == 0) {
                assert(meta_point.size() > s_block);
                meta_point[s_block] = offset;
                metap_cont.push_back(offset);
                assert(meta_rank.size() > s_block);
                meta_rank[s_block++] = sum_;
                metar_cont.push_back(sum_);
            }
            
            uint64_t elem = bv.get_int(block * block_width, block_width);
            auto enc = coder.encode(elem);
            /*if (block == 269) {
                std::cerr << block << " block p_off: " << meta_point[s_block - 1] << "\n"
                          << std::bitset<64>(elem) << ", " << elem << "\n"
                          << " -> " << enc.first << ", " << enc.second << "\n"
                          << " actual offset = " << offset << "\n"
                          << " from " << block * block_width << std::endl;
            }*/
            // std::cerr << std::bitset<64>(elem) << " -> " << enc.first << ", "
            // << enc.second << std::endl;
            sum_ += enc.first;
            assert(data_typ.size() > block);
            data_typ[block++] = enc.first;
            datat_cont.push_back(enc.first);
            data_val.set_int(offset, enc.second, widths[enc.first]);
            datav_cont.push_back(enc.second);
            // std::cerr << "                -> " << widths[enc.first] <<
            // std::endl;
            offset += widths[enc.first];
        }
        /*for (uint64_t i = 0; i < metap_cont.size(); ++i) {
            if (metap_cont[i] != meta_point[i]) {
                std::cerr << "s_block " << i << " offset pointer missmatch. Should be " << metap_cont[i] << " and is " << meta_point[i] << std::endl;
                exit(1);
            }
            if (metar_cont[i] != meta_rank[i]) {
                std::cerr << "s_block " << i << " partial rank missmatch. Should be " << metar_cont[i] << " and is " << meta_rank[i] << std::endl;
                exit(1);
            }
        }
        uint64_t off_cont = 0;
        for (uint64_t i = 0; i < datat_cont.size(); ++i) {
            if (i % k == 0) {
                if (meta_point[i / k] != off_cont) {
                    std::cerr << "s_block " << (i / k) << " offset pointer missmarch. calculated: " << off_cont << ", stored: " << meta_point[i / k] << std::endl;
                    exit(1);
                }
            }
            if (datat_cont[i] != data_typ[i]) {
                std::cerr << "block " << i << " data type missmatch. Should be " << datat_cont[i] << " and is " << data_typ[i] << std::endl;
                exit(1);
            }
            uint16_t w = widths[datat_cont[i]];
            uint64_t f_read = data_val.get_int(off_cont, w);
            if (f_read != datav_cont[i]) {
                std::cerr << "block " << i << " data value missmatch: Width = " << w << ", Class = " << datat_cont[i] << "\n"
                          << std::bitset<64>(f_read) << ", " << f_read << "\n"
                          << std::bitset<64>(datav_cont[i]) << ", " << datav_cont[i] << std::endl;
                exit(1);
            }
            auto enc = coder.encode(bv.get_int(i * 63, 63));
            if (datat_cont[i] != enc.first || datav_cont[i] != enc.second) {
                std::cerr << "WTF encoding pass difference... (" << enc.first << ", " << enc.second << ") != (" << datat_cont[i] << ", " << datav_cont[i] << ")" << std::endl;
            }
            off_cont += w;
        }
        std::cerr << "data check OK" << std::endl;*/
    }

    bool access(uint64_t i) const {
        uint64_t block = i / block_width;
        uint64_t s_block = block / k;
        uint16_t b_offset = block % k;
        i %= block_width;
        uint64_t offset = meta_point[s_block];
        /*if (block == 269) {
            std::cerr << block << " block s_offset: " << offset << std::endl;
        }*/
        for (uint64_t j = 0; j < k; ++j) {
            if (b_offset == j) {
                break;
            }
            auto bk = data_typ[s_block * k + j];
            offset += widths[bk];
        }
        /*if (block == 269) {
            std::cerr << block << " block actual offset: " << offset << std::endl;
        }*/
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