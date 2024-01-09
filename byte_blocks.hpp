
#include <cstdint>
#include <array>

template <uint16_t block_size>
class rrr_decode {
   private:
    static const constexpr uint16_t MAX_BH = (1 << 8) - 1;

    struct f_stop {
        uint64_t f;
        uint8_t byte;
        uint8_t pop;
    };

    static const constexpr std::array<std::array<f_stop, b>, b / 8> f_stops;
    

    inline uint32_t find(uint16_t bs, uint16_t k, uint64_t f) {
        uint16_t a = 0;
        uint16_t b = MAX_BH;
        while (a < b) {
            uint16_t m = (a + b + 1) / 2;
            if (f_stops[bs][k][m].f > f) {
                b = m - 1;
            } else {
                a = m;
            }
        }
        return a;
    }

    template <uint16_t bs>
    bool access(uint16_t k, uint64_t f, uint16_t i) {
        uint32_t boff = find(bs, k, f);
        auto block = f_stops[bs][k][boff];
        if (bs - i < 8) {
            return (block.byte >> (bs - i)) & 1;
        }
        access<bs - 8>(k - block.pop, f - block.f, i - 8);
    }

   public:
    bool access(uint16_t k, uint64_t f, uint16_t i) {
        return access<block_size>(k, f, i);
    }
};