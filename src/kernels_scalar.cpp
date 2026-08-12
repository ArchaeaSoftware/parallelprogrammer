#include <cstring>

#include "cbfp/limb_column.hpp"
#include "kernels.hpp"

namespace cbfp {
namespace kernels {
namespace {

struct Split {
    std::uint64_t mantissa;  // odd, or 0 for +/-0
    long long exponent;      // true-ulp exponent
    long long top;           // one past the highest bit position reached
    bool negative;
    bool nonfinite;
};

// The IEEE-754 field extraction both passes share, kept in one place so the
// scan and the accumulate cannot disagree about what a double means.
inline Split split(double v)
{
    std::uint64_t bits;
    std::memcpy(&bits, &v, sizeof bits);

    const std::uint64_t biased = (bits >> 52) & 0x7FF;
    const std::uint64_t frac = bits & ((std::uint64_t{1} << 52) - 1);

    std::uint64_t m = frac;
    if (biased != 0) m |= std::uint64_t{1} << 52;
    long long e = static_cast<long long>(biased < 1 ? 1 : biased) - 1075;

    Split s;
    s.negative = (bits >> 63) != 0;
    s.nonfinite = (biased == 0x7FF);
    // Normalizing to odd raises the exponent by the trailing zero count and
    // lowers the significand's width by the same amount, so `top` is
    // unaffected and the scan never needs the count.
    s.top = e + (m == 0 ? 0 : 64 - __builtin_clzll(m));
    if (m != 0) {
        const int tz = __builtin_ctzll(m);
        m >>= tz;
        e += tz;
    }
    s.mantissa = m;
    s.exponent = e;
    return s;
}

}  // namespace

Scan scan_column_scalar(const double* values, std::size_t rows,
                        long long floor_exponent)
{
    // First pass reads only the exponent field. The raw exponent is a lower
    // bound on the true ulp, so if it clears the floor no rescale is due and
    // the significand never has to be examined.
    Scan out{0, 0, false, false};
    long long min_raw = 0;
    std::uint64_t max_abs = 0;

    for (std::size_t i = 0; i < rows; ++i) {
        std::uint64_t bits;
        std::memcpy(&bits, &values[i], sizeof bits);
        const std::uint64_t abs_bits = bits & 0x7FFFFFFFFFFFFFFFull;
        const std::uint64_t biased = abs_bits >> 52;
        if (biased == 0x7FF) {
            out.nonfinite = true;
            return out;
        }
        if (abs_bits == 0) continue;
        const long long e =
            static_cast<long long>(biased < 1 ? 1 : biased) - 1075;
        if (!out.any) {
            min_raw = e;
            out.any = true;
        } else if (e < min_raw) {
            min_raw = e;
        }
        // With the sign cleared, the bit pattern orders by magnitude, so the
        // largest pattern is the largest value.
        if (abs_bits > max_abs) max_abs = abs_bits;
    }
    if (!out.any) return out;

    const std::uint64_t max_biased = max_abs >> 52;
    out.max_top = max_biased != 0 ? static_cast<long long>(max_biased) - 1022
                                  : -1074 + (64 - __builtin_clzll(max_abs));

    if (min_raw >= floor_exponent) {
        out.min_exponent = min_raw;
        return out;
    }

    // Only now, when the exponent could actually drop, is the exact true-ulp
    // minimum worth computing.
    bool first = true;
    for (std::size_t i = 0; i < rows; ++i) {
        const Split s = split(values[i]);
        if (s.mantissa == 0) continue;
        if (first || s.exponent < out.min_exponent) {
            out.min_exponent = s.exponent;
            first = false;
        }
    }
    return out;
}

void accumulate_one(std::uint64_t* const* limbs, std::size_t nlimbs,
                    std::size_t row, std::uint64_t mantissa, std::size_t shift,
                    bool negative)
{
    if (mantissa == 0) return;

    const std::size_t off = shift / 64;
    const unsigned bit = shift % 64;
    const std::uint64_t lo = mantissa << bit;
    const std::uint64_t hi = bit != 0 ? (mantissa >> (64 - bit)) : 0;

    if (negative) {
        std::uint64_t borrow = 0;
        for (std::size_t p = off; p < nlimbs; ++p) {
            const std::uint64_t a = (p == off) ? lo : ((p == off + 1) ? hi : 0);
            const std::uint64_t x = limbs[p][row];
            const std::uint64_t d = x - a;
            const std::uint64_t b1 = (x < a) ? 1 : 0;
            const std::uint64_t d2 = d - borrow;
            const std::uint64_t b2 = (d < borrow) ? 1 : 0;
            limbs[p][row] = d2;
            borrow = b1 | b2;
            // Past the addend with no borrow left, nothing further changes.
            if (borrow == 0 && p >= off + 1) break;
        }
    } else {
        std::uint64_t carry = 0;
        for (std::size_t p = off; p < nlimbs; ++p) {
            const std::uint64_t a = (p == off) ? lo : ((p == off + 1) ? hi : 0);
            const std::uint64_t x = limbs[p][row];
            const std::uint64_t s = x + a;
            const std::uint64_t c1 = (s < x) ? 1 : 0;
            const std::uint64_t s2 = s + carry;
            const std::uint64_t c2 = (s2 < s) ? 1 : 0;
            limbs[p][row] = s2;
            carry = c1 | c2;
            if (carry == 0 && p >= off + 1) break;
        }
    }
}

void accumulate_scalar(std::uint64_t* const* limbs, std::size_t nlimbs,
                       const double* values, std::size_t rows,
                       std::int32_t column_exponent, std::size_t first_limb)
{
    (void)first_limb;  // the scalar path starts from each row's own limb
    for (std::size_t i = 0; i < rows; ++i) {
        const Split s = split(values[i]);
        if (s.mantissa == 0) continue;
        accumulate_one(limbs, nlimbs, i, s.mantissa,
                       static_cast<std::size_t>(s.exponent - column_exponent),
                       s.negative);
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
