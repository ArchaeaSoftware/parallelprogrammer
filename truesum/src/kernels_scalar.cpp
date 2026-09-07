#include <cstring>

#include "truesum/limb_column.hpp"
#include "kernels.hpp"

namespace truesum {
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
// survey and the accumulate cannot disagree about what a double means.
inline Split
split(double v)
{
    std::uint64_t bits;
    std::memcpy(&bits, &v, sizeof bits);

    const std::uint64_t biased = (bits >> 52) & 0x7FF;
    const std::uint64_t frac = bits & ((std::uint64_t{1} << 52) - 1);

    std::uint64_t m = frac;
    if (0 != biased) m |= std::uint64_t{1} << 52;
    long long e = static_cast<long long>(biased < 1 ? 1 : biased) - 1075;

    Split s;
    s.negative = (bits >> 63) != 0;
    s.nonfinite = (biased == 0x7FF);
    // Normalizing to odd raises the exponent by the trailing zero count and
    // lowers the significand's width by the same amount, so `top` is
    // unaffected and the survey never needs the count.
    s.top = e + (0 == m ? 0 : 64 - __builtin_clzll(m));
    if (0 != m) {
        const int tz = __builtin_ctzll(m);
        m >>= tz;
        e += tz;
    }
    s.mantissa = m;
    s.exponent = e;
    return s;
}

}  // namespace

Survey
survey_column_scalar(const double *values, std::size_t rows,
                     long long floor_exponent)
{
    // First pass reads only the exponent field. The raw exponent is a lower
    // bound on the true ulp, so if it clears the floor no rescale is due and
    // the significand never has to be examined.
    Survey out{0, 0, false, false};
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
        if (0 == abs_bits) continue;
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
    out.max_top = 0 != max_biased ? static_cast<int>(max_biased) - 1022
                                  : -1074 + (64 - __builtin_clzll(max_abs));

    if (min_raw >= floor_exponent) {
        out.min_exponent = static_cast<int>(min_raw);
        return out;
    }

    // Only now, when the exponent could actually drop, is the exact true-ulp
    // minimum worth computing.
    bool first = true;
    for (std::size_t i = 0; i < rows; ++i) {
        const Split s = split(values[i]);
        if (0 == s.mantissa) continue;
        if (first || s.exponent < out.min_exponent) {
            out.min_exponent = static_cast<int>(s.exponent);
            first = false;
        }
    }
    return out;
}

// One limb of the carry chain, both directions. `carry` is a borrow when
// subtracting. Written once here so the loop below and its peeled top limb
// cannot drift apart.
inline std::uint64_t
add_limb(std::uint64_t x, std::uint64_t a, std::uint64_t &carry)
{
    const std::uint64_t s = x + a;
    const std::uint64_t c1 = (s < x) ? 1 : 0;
    const std::uint64_t s2 = s + carry;
    const std::uint64_t c2 = (s2 < s) ? 1 : 0;
    carry = c1 | c2;
    return s2;
}

inline std::uint64_t
sub_limb(std::uint64_t x, std::uint64_t a, std::uint64_t &borrow)
{
    const std::uint64_t d = x - a;
    const std::uint64_t b1 = (x < a) ? 1 : 0;
    const std::uint64_t d2 = d - borrow;
    const std::uint64_t b2 = (d < borrow) ? 1 : 0;
    borrow = b1 | b2;
    return d2;
}

// One row's limbs are loaded once, every batch's addend applied to them in
// registers, and the result written once. The inner loop is the same carry
// chain accumulate_one runs, with `v` standing in for the limb arrays.
void
accumulate_fold_scalar(std::uint64_t *const *limbs, std::size_t nlimbs,
                       const double *const *columns, std::size_t count,
                       std::size_t rows, std::int32_t column_exponent,
                       unsigned *flags)
{
    unsigned bad = 0;
    std::uint64_t overflow = 0;
    std::uint64_t v[kMaxFoldLimbs];

    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t k = 0; k < nlimbs; ++k) v[k] = limbs[k][i];

        for (std::size_t b = 0; b < count; ++b) {
            const Split s = split(columns[b][i]);
            if (s.nonfinite) bad |= kBadNonFinite;
            if (0 == s.mantissa) continue;
            const long long shift =
                s.exponent - static_cast<long long>(column_exponent);
            if (shift < 0) {
                bad |= kBadExponent;
                continue;
            }
            // An addend must stop below the sign bit; one that reaches it
            // has no representation at this width, and would also defeat
            // the sign test at the top limb below.
            if (s.top - static_cast<long long>(column_exponent) >=
                static_cast<long long>(64 * nlimbs)) {
                bad |= kBadWidth;
                continue;
            }
            const std::size_t off = static_cast<std::size_t>(shift) / 64;
            const unsigned bit = static_cast<unsigned>(shift) % 64;
            const std::uint64_t lo =
                0 == bit ? s.mantissa : s.mantissa << bit;
            const std::uint64_t hi =
                0 == bit ? 0 : s.mantissa >> (64 - bit);

            // Sign test after the chain, on the last limb it touched, and
            // only when that was the top one; see accumulate_one for why.
            std::uint64_t x = 0, r = 0, carry = 0;
            std::size_t p = off;
            if (s.negative) {
                for (; p < nlimbs; ++p) {
                    const std::uint64_t a =
                        (p == off) ? lo : ((p == off + 1) ? hi : 0);
                    x = v[p];
                    r = sub_limb(x, a, carry);
                    v[p] = r;
                    if (0 == carry && p >= off + 1) break;
                }
                if (p + 1 >= nlimbs) overflow |= (x & ~r) >> 63;
            } else {
                for (; p < nlimbs; ++p) {
                    const std::uint64_t a =
                        (p == off) ? lo : ((p == off + 1) ? hi : 0);
                    x = v[p];
                    r = add_limb(x, a, carry);
                    v[p] = r;
                    if (0 == carry && p >= off + 1) break;
                }
                if (p + 1 >= nlimbs) overflow |= (~x & r) >> 63;
            }
        }

        for (std::size_t k = 0; k < nlimbs; ++k) limbs[k][i] = v[k];
    }
    if (0 != overflow) bad |= kBadOverflow;
    *flags |= bad;
}

// Returns true if the sum overflowed its width, the addend having fit. The
// caller has already checked that it fits; see accumulate_scalar.
//
// The sign test is applied after the chain, to the last limb it touched, and
// only when that limb was the top one. Testing inside the loop body cost the
// two-limb case 45% when measured, and peeling the top limb out of the loop
// still cost the three-limb case 27%; hoisting the last iteration's operands
// out of the loop leaves the loop as it was.
bool
accumulate_one(std::uint64_t *const *limbs, std::size_t nlimbs, std::size_t row,
               std::uint64_t mantissa, std::size_t shift, bool negative)
{
    if (0 == mantissa) return false;

    const std::size_t off = shift / 64;
    const unsigned bit = shift % 64;
    const std::uint64_t lo = mantissa << bit;
    const std::uint64_t hi = 0 != bit ? (mantissa >> (64 - bit)) : 0;

    std::uint64_t x = 0, r = 0;  // the last limb's old and new value
    std::size_t p = off;
    if (negative) {
        std::uint64_t borrow = 0;
        for (; p < nlimbs; ++p) {
            const std::uint64_t a = (p == off) ? lo : ((p == off + 1) ? hi : 0);
            x = limbs[p][row];
            r = sub_limb(x, a, borrow);
            limbs[p][row] = r;
            // Past the addend with no borrow left, nothing further changes.
            if (0 == borrow && p >= off + 1) break;
        }
        // A negative value that lost its sign bit overflowed.
        return p + 1 >= nlimbs && 0 != ((x & ~r) >> 63);
    }
    std::uint64_t carry = 0;
    for (; p < nlimbs; ++p) {
        const std::uint64_t a = (p == off) ? lo : ((p == off + 1) ? hi : 0);
        x = limbs[p][row];
        r = add_limb(x, a, carry);
        limbs[p][row] = r;
        if (0 == carry && p >= off + 1) break;
    }
    // A non-negative value that gained a sign bit overflowed.
    return p + 1 >= nlimbs && 0 != ((~x & r) >> 63);
}

void
accumulate_scalar(std::uint64_t *const *limbs, std::size_t nlimbs,
                  const double *values, std::size_t rows,
                  std::int32_t column_exponent, std::size_t first_limb,
                  unsigned *flags)
{
    (void)first_limb;  // the scalar path starts from each row's own limb
    unsigned bad = 0;
    for (std::size_t i = 0; i < rows; ++i) {
        const Split s = split(values[i]);
        if (s.nonfinite) bad |= kBadNonFinite;
        if (0 == s.mantissa) continue;
        const long long shift =
            s.exponent - static_cast<long long>(column_exponent);
        if (shift < 0) {
            bad |= kBadExponent;
            continue;
        }
        if (s.top - static_cast<long long>(column_exponent) >=
            static_cast<long long>(64 * nlimbs)) {
            bad |= kBadWidth;
            continue;
        }
        if (accumulate_one(limbs, nlimbs, i, s.mantissa,
                           static_cast<std::size_t>(shift), s.negative)) {
            bad |= kBadOverflow;
        }
    }
    *flags |= bad;
}

void
shift_left(std::uint64_t *const *dst, std::size_t ndst,
           std::uint64_t *const *src, std::size_t nsrc, std::size_t rows,
           unsigned shift)
{
    const std::size_t word = shift / 64;
    const unsigned bit = shift % 64;
    const std::size_t n = padded_rows(rows);

    for (std::size_t k = ndst; k-- > 0;) {
        std::uint64_t *out = dst[k];
        if (k < word) {
            for (std::size_t i = 0; i < n; ++i) out[i] = 0;
            continue;
        }
        const std::size_t a = k - word;  // contributes its low part
        const std::size_t c = a - 1;     // contributes its high part

        // Above the source, every limb reads as the sign fill.
        if (a >= nsrc) {
            const std::uint64_t *top = src[nsrc - 1];
            if (0 == bit || c >= nsrc) {
                for (std::size_t i = 0; i < n; ++i) {
                    out[i] = static_cast<std::uint64_t>(
                        static_cast<std::int64_t>(top[i]) >> 63);
                }
            } else {
                const std::uint64_t *low = src[c];
                for (std::size_t i = 0; i < n; ++i) {
                    const std::uint64_t fill = static_cast<std::uint64_t>(
                        static_cast<std::int64_t>(top[i]) >> 63);
                    out[i] = (fill << bit) | (low[i] >> (64 - bit));
                }
            }
            continue;
        }

        const std::uint64_t *high = src[a];
        if (0 == bit) {
            for (std::size_t i = 0; i < n; ++i) out[i] = high[i];
        } else if (0 == a) {
            for (std::size_t i = 0; i < n; ++i) out[i] = high[i] << bit;
        } else {
            const std::uint64_t *low = src[c];
            for (std::size_t i = 0; i < n; ++i) {
                out[i] = (high[i] << bit) | (low[i] >> (64 - bit));
            }
        }
    }
}

void
sign_fill(std::uint64_t *dst, const std::uint64_t *top, std::size_t rows)
{
    const std::size_t n = padded_rows(rows);
    for (std::size_t i = 0; i < n; ++i) {
        dst[i] =
            static_cast<std::uint64_t>(static_cast<std::int64_t>(top[i]) >> 63);
    }
}

}  // namespace kernels
}  // namespace truesum
