#include "cbfp/limb_column.hpp"
#include "kernels.hpp"

namespace cbfp {
namespace kernels {

void accumulate_scalar(std::uint64_t* const* limbs, std::size_t nlimbs,
                       const Batch& b)
{
    for (std::size_t i = 0; i < b.rows; ++i) {
        const std::uint64_t m = b.mantissa[i];
        if (m == 0) continue;

        const std::size_t k =
            static_cast<std::size_t>(b.exponent[i] - b.column_exponent);
        const std::size_t off = k / 64;
        const unsigned bit = k % 64;
        const std::uint64_t lo = m << bit;
        const std::uint64_t hi = bit != 0 ? (m >> (64 - bit)) : 0;

        if (b.negative[i] != 0) {
            std::uint64_t borrow = 0;
            for (std::size_t p = off; p < nlimbs; ++p) {
                const std::uint64_t a =
                    (p == off) ? lo : ((p == off + 1) ? hi : 0);
                const std::uint64_t x = limbs[p][i];
                const std::uint64_t d = x - a;
                const std::uint64_t b1 = (x < a) ? 1 : 0;
                const std::uint64_t d2 = d - borrow;
                const std::uint64_t b2 = (d < borrow) ? 1 : 0;
                limbs[p][i] = d2;
                borrow = b1 | b2;
                // Past the addend with no borrow left, nothing further changes.
                if (borrow == 0 && p >= off + 1) break;
            }
        } else {
            std::uint64_t carry = 0;
            for (std::size_t p = off; p < nlimbs; ++p) {
                const std::uint64_t a =
                    (p == off) ? lo : ((p == off + 1) ? hi : 0);
                const std::uint64_t x = limbs[p][i];
                const std::uint64_t s = x + a;
                const std::uint64_t c1 = (s < x) ? 1 : 0;
                const std::uint64_t s2 = s + carry;
                const std::uint64_t c2 = (s2 < s) ? 1 : 0;
                limbs[p][i] = s2;
                carry = c1 | c2;
                if (carry == 0 && p >= off + 1) break;
            }
        }
    }
}

void shift_left(std::uint64_t* const* dst, std::size_t ndst,
                std::uint64_t* const* src, std::size_t nsrc, std::size_t rows,
                unsigned shift)
{
    const std::size_t word = shift / 64;
    const unsigned bit = shift % 64;
    const std::size_t n = padded_rows(rows);

    for (std::size_t k = ndst; k-- > 0;) {
        std::uint64_t* out = dst[k];
        if (k < word) {
            for (std::size_t i = 0; i < n; ++i) out[i] = 0;
            continue;
        }
        const std::size_t a = k - word;  // contributes its low part
        const std::size_t c = a - 1;     // contributes its high part

        // Above the source, every limb reads as the sign fill.
        if (a >= nsrc) {
            const std::uint64_t* top = src[nsrc - 1];
            if (bit == 0 || c >= nsrc) {
                for (std::size_t i = 0; i < n; ++i) {
                    out[i] = static_cast<std::uint64_t>(
                        static_cast<std::int64_t>(top[i]) >> 63);
                }
            } else {
                const std::uint64_t* low = src[c];
                for (std::size_t i = 0; i < n; ++i) {
                    const std::uint64_t fill = static_cast<std::uint64_t>(
                        static_cast<std::int64_t>(top[i]) >> 63);
                    out[i] = (fill << bit) | (low[i] >> (64 - bit));
                }
            }
            continue;
        }

        const std::uint64_t* high = src[a];
        if (bit == 0) {
            for (std::size_t i = 0; i < n; ++i) out[i] = high[i];
        } else if (a == 0) {
            for (std::size_t i = 0; i < n; ++i) out[i] = high[i] << bit;
        } else {
            const std::uint64_t* low = src[c];
            for (std::size_t i = 0; i < n; ++i) {
                out[i] = (high[i] << bit) | (low[i] >> (64 - bit));
            }
        }
    }
}

void sign_fill(std::uint64_t* dst, const std::uint64_t* top, std::size_t rows)
{
    const std::size_t n = padded_rows(rows);
    for (std::size_t i = 0; i < n; ++i) {
        dst[i] =
            static_cast<std::uint64_t>(static_cast<std::int64_t>(top[i]) >> 63);
    }
}

}  // namespace kernels
}  // namespace cbfp
