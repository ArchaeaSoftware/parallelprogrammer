// Low-level arbitrary-precision limb arithmetic.
//
// A "value" is a little-endian array of 64-bit limbs interpreted as a two's
// complement signed integer. All routines here are magnitude/width agnostic:
// the caller owns allocation and is responsible for providing enough limbs.
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace cbfp {
namespace limbs {

using limb_t = std::uint64_t;
inline constexpr std::size_t kLimbBits = 64;

// True if the two's complement value is negative (top bit of the top limb).
bool
is_negative(const limb_t* p, std::size_t n);

bool
is_zero(const limb_t* p, std::size_t n);

// dst = src, sign-extended to dn limbs. Requires dn >= sn.
void
widen(const limb_t* src, std::size_t sn, limb_t* dst, std::size_t dn);

// dst = src << s, sign-aware. Requires dn limbs to be enough to hold the
// result without losing significant bits; the caller sizes dn.
void
shift_left(const limb_t* src, std::size_t sn, unsigned s, limb_t* dst,
           std::size_t dn);

// p += (hi:lo) << (64*off), or p -= that when subtract is true.
// (hi:lo) is an unsigned 128-bit magnitude. Returns true on signed overflow,
// in which case p holds the wrapped result and the caller must widen + retry.
bool
add_shifted(limb_t* p, std::size_t n, limb_t lo, limb_t hi, std::size_t off,
            bool subtract);

// dst = -src (two's complement negation), n limbs.
void
negate_into(const limb_t* src, std::size_t n, limb_t* dst);

// Number of significant bits of a NON-NEGATIVE value (0 for zero).
std::size_t
bit_length(const limb_t* p, std::size_t n);

bool
get_bit(const limb_t* p, std::size_t n, std::size_t i);

// True if any bit strictly below index i is set.
bool
any_bits_below(const limb_t* p, std::size_t n, std::size_t i);

// Low 64 bits of (p >> shift), for a non-negative p.
std::uint64_t
extract_u64(const limb_t* p, std::size_t n, std::size_t shift);

// v = v * m (unsigned, grows v as needed).
void
mul_small(std::vector<limb_t>& v, limb_t m);

// Exact decimal digits of a non-negative magnitude. Returns "0" for zero.
std::string
magnitude_to_decimal(std::vector<limb_t> mag);

}  // namespace limbs
}  // namespace cbfp
