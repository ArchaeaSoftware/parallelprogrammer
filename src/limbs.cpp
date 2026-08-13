#include "cbfp/limbs.hpp"

#include <cassert>

namespace cbfp {
namespace limbs {
namespace {

// Sign-extended limb read: indices at or past the top return the sign fill.
inline limb_t
at(const limb_t *p, std::size_t n, std::size_t i, limb_t fill)
{
    return i < n ? p[i] : fill;
}

// Drops leading zero limbs from a magnitude.
void
trim(std::vector<limb_t> &v)
{
    std::size_t n = v.size();
    for (; n > 0 && 0 == v[n - 1]; --n) {
    }
    v.resize(n);
}

}  // namespace

bool
is_negative(const limb_t *p, std::size_t n)
{
    return n != 0 && (p[n - 1] >> 63) != 0;
}

bool
is_zero(const limb_t *p, std::size_t n)
{
    for (std::size_t i = 0; i < n; ++i) {
        if (0 != p[i]) return false;
    }
    return true;
}

void
widen(const limb_t *src, std::size_t sn, limb_t *dst, std::size_t dn)
{
    assert(dn >= sn);
    const limb_t fill = is_negative(src, sn) ? ~limb_t{0} : limb_t{0};
    for (std::size_t i = 0; i < sn; ++i) dst[i] = src[i];
    for (std::size_t i = sn; i < dn; ++i) dst[i] = fill;
}

void
shift_left(const limb_t *src, std::size_t sn, unsigned s, limb_t *dst,
           std::size_t dn)
{
    const limb_t fill = is_negative(src, sn) ? ~limb_t{0} : limb_t{0};
    const std::size_t word = s / kLimbBits;
    const unsigned bit = s % kLimbBits;

    for (std::size_t i = dn; i-- > 0;) {
        if (i < word) {
            dst[i] = 0;
            continue;
        }
        const std::size_t k = i - word;
        limb_t v = at(src, sn, k, fill) << bit;
        if (0 != bit && k > 0) {
            v |= at(src, sn, k - 1, fill) >> (kLimbBits - bit);
        }
        dst[i] = v;
    }
}

bool
add_shifted(limb_t *p, std::size_t n, limb_t lo, limb_t hi, std::size_t off,
            bool subtract)
{
    const bool was_negative = is_negative(p, n);
    const limb_t addend[2] = {lo, hi};

    if (subtract) {
        limb_t borrow = 0;
        for (std::size_t i = 0; i < n; ++i) {
            const limb_t a = (i >= off && i < off + 2) ? addend[i - off] : 0;
            const __uint128_t d = static_cast<__uint128_t>(p[i]) - a - borrow;
            p[i] = static_cast<limb_t>(d);
            borrow = static_cast<limb_t>(d >> 127) & 1;
        }
        // Subtracting a non-negative value from a negative one must stay
        // negative.
        return was_negative && !is_negative(p, n);
    }

    limb_t carry = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const limb_t a = (i >= off && i < off + 2) ? addend[i - off] : 0;
        const __uint128_t s = static_cast<__uint128_t>(p[i]) + a + carry;
        p[i] = static_cast<limb_t>(s);
        carry = static_cast<limb_t>(s >> 64);
    }
    // Adding a non-negative value to a non-negative one must stay
    // non-negative.
    return !was_negative && is_negative(p, n);
}

void
negate_into(const limb_t *src, std::size_t n, limb_t *dst)
{
    limb_t carry = 1;
    for (std::size_t i = 0; i < n; ++i) {
        const __uint128_t s = static_cast<__uint128_t>(~src[i]) + carry;
        dst[i] = static_cast<limb_t>(s);
        carry = static_cast<limb_t>(s >> 64);
    }
}

std::size_t
bit_length(const limb_t *p, std::size_t n)
{
    for (std::size_t i = n; i-- > 0;) {
        if (0 != p[i]) {
            return i * kLimbBits + (kLimbBits - __builtin_clzll(p[i]));
        }
    }
    return 0;
}

bool
get_bit(const limb_t *p, std::size_t n, std::size_t i)
{
    const std::size_t w = i / kLimbBits;
    if (w >= n) return false;
    return (p[w] >> (i % kLimbBits)) & 1;
}

bool
any_bits_below(const limb_t *p, std::size_t n, std::size_t i)
{
    const std::size_t w = i / kLimbBits;
    const unsigned b = i % kLimbBits;
    for (std::size_t k = 0; k < w && k < n; ++k) {
        if (0 != p[k]) return true;
    }
    if (0 != b && w < n) {
        if (0 != (p[w] & ((limb_t{1} << b) - 1))) return true;
    }
    return false;
}

std::uint64_t
extract_u64(const limb_t *p, std::size_t n, std::size_t shift)
{
    const std::size_t w = shift / kLimbBits;
    const unsigned b = shift % kLimbBits;
    limb_t v = at(p, n, w, 0) >> b;
    if (0 != b) {
        v |= at(p, n, w + 1, 0) << (kLimbBits - b);
    }
    return v;
}

void
mul_small(std::vector<limb_t> &v, limb_t m)
{
    limb_t carry = 0;
    for (std::size_t i = 0; i < v.size(); ++i) {
        const __uint128_t prod = static_cast<__uint128_t>(v[i]) * m + carry;
        v[i] = static_cast<limb_t>(prod);
        carry = static_cast<limb_t>(prod >> 64);
    }
    if (0 != carry) v.push_back(carry);
}

std::string
magnitude_to_decimal(std::vector<limb_t> mag)
{
    trim(mag);
    if (mag.empty()) return "0";

    constexpr limb_t kChunk = 10000000000000000000ull;  // 10^19
    constexpr std::size_t kChunkDigits = 19;

    // Each division by 10^19 > 2^63 removes at least 63 bits, which bounds the
    // number of chunks the value can produce.
    const std::size_t bits = bit_length(mag.data(), mag.size());
    const std::size_t max_chunks = bits / 63 + 1;

    std::vector<std::string> chunks;
    chunks.reserve(max_chunks);
    for (std::size_t c = 0; c < max_chunks && !mag.empty(); ++c) {
        limb_t rem = 0;
        for (std::size_t i = mag.size(); i-- > 0;) {
            const __uint128_t cur =
                (static_cast<__uint128_t>(rem) << 64) | mag[i];
            mag[i] = static_cast<limb_t>(cur / kChunk);
            rem = static_cast<limb_t>(cur % kChunk);
        }
        trim(mag);
        chunks.push_back(std::to_string(rem));
    }
    assert(mag.empty());

    std::string out = chunks.back();
    for (std::size_t i = chunks.size() - 1; i-- > 0;) {
        out.append(kChunkDigits - chunks[i].size(), '0');
        out.append(chunks[i]);
    }
    return out;
}

}  // namespace limbs
}  // namespace cbfp
