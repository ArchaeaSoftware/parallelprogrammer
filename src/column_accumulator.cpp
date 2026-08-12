#include "cbfp/column_accumulator.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace cbfp {

using limbs::kLimbBits;
using limbs::limb_t;

namespace {

std::size_t limbs_for_bits(std::size_t bits)
{
    return (bits + kLimbBits - 1) / kLimbBits;
}

// Number of trailing zero bits of a non-negative magnitude (0 if zero).
std::size_t trailing_zeros(const std::vector<limb_t>& v)
{
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (v[i] != 0) return i * kLimbBits + __builtin_ctzll(v[i]);
    }
    return 0;
}

std::size_t ceil_log2(std::size_t n)
{
    if (n <= 1) return 0;
    return kLimbBits - __builtin_clzll(n - 1);
}

// 5^k for k <= 27, which is the largest power of five that fits in a limb.
limb_t pow5(std::size_t k)
{
    assert(k <= 27);
    limb_t p = 1;
    for (std::size_t n = 0; n < k; ++n) p *= 5;
    return p;
}

}  // namespace

DoubleParts decompose(double v)
{
    DoubleParts p{0, false, 0};
    if (v == 0.0) return p;  // also catches -0.0

    p.negative = std::signbit(v);
    int e = 0;
    const double f =
        std::frexp(std::fabs(v), &e);  // |v| = f * 2^e, f in [.5,1)
    std::uint64_t m = static_cast<std::uint64_t>(std::ldexp(f, 53));  // exact
    e -= 53;

    // Normalizing the mantissa to odd makes the exponent the value's true ulp,
    // which keeps columns of coarse values (integers, say) narrow.
    const int tz = __builtin_ctzll(m);
    m >>= tz;
    e += tz;

    p.mantissa = m;
    p.exponent = e;
    return p;
}

ColumnBlockMatrix::ColumnBlockMatrix(std::size_t rows, std::size_t cols)
    : rows_(rows), cols_(cols), cols_state_(cols)
{
    for (auto& c : cols_state_) c.data.assign(rows_, 0);
}

void ColumnBlockMatrix::check_index(std::size_t i, std::size_t j) const
{
    if (i >= rows_ || j >= cols_) {
        throw std::out_of_range("cbfp: matrix index out of range");
    }
}

void ColumnBlockMatrix::ensure_limbs(Column& c, std::size_t limbs_needed)
{
    if (c.limbs >= limbs_needed) return;

    std::vector<limb_t> wider(rows_ * limbs_needed);
    for (std::size_t i = 0; i < rows_; ++i) {
        limbs::widen(c.data.data() + i * c.limbs, c.limbs,
                     wider.data() + i * limbs_needed, limbs_needed);
    }
    c.data.swap(wider);
    c.limbs = limbs_needed;
}

void ColumnBlockMatrix::rescale(Column& c, int new_exponent)
{
    if (new_exponent >= c.exponent) return;
    const unsigned shift = static_cast<unsigned>(
        static_cast<long long>(c.exponent) - new_exponent);

    // An all-zero column carries no information, so its scale can move freely
    // without widening anything.
    if (limbs::is_zero(c.data.data(), c.data.size())) {
        c.exponent = new_exponent;
        return;
    }

    const std::size_t wide = limbs_for_bits(c.limbs * kLimbBits + shift);
    std::vector<limb_t> shifted(rows_ * wide);
    for (std::size_t i = 0; i < rows_; ++i) {
        limbs::shift_left(c.data.data() + i * c.limbs, c.limbs, shift,
                          shifted.data() + i * wide, wide);
    }
    c.data.swap(shifted);
    c.limbs = wide;
    c.exponent = new_exponent;
}

void ColumnBlockMatrix::accumulate(std::size_t i, std::size_t j, double v,
                                   bool negate, int log2_scale)
{
    check_index(i, j);
    if (!std::isfinite(v)) {
        throw std::domain_error("cbfp: cannot accumulate a non-finite value");
    }

    const DoubleParts p = decompose(v);
    if (p.mantissa == 0) return;

    const bool subtract = (p.negative != negate);
    const long long e64 = static_cast<long long>(p.exponent) + log2_scale;
    // A double's own exponent is bounded by +/-1074, so only an extreme
    // log2_scale can push this out of range; refuse rather than wrap.
    if (e64 < -(1 << 24) || e64 > (1 << 24)) {
        throw std::domain_error("cbfp: log2_scale puts the value out of range");
    }
    const int e = static_cast<int>(e64);

    Column& c = cols_state_[j];
    if (!c.initialized) {
        c.exponent = e;
        c.initialized = true;
    }
    if (e < c.exponent) rescale(c, e);

    const std::size_t k = static_cast<std::size_t>(e - c.exponent);
    const std::size_t off = k / kLimbBits;
    const unsigned bit = k % kLimbBits;

    // k + 53 magnitude bits + 1 sign bit must fit before the add is attempted;
    // overflow past that is handled below.
    ensure_limbs(c, limbs_for_bits(k + 54));

    const limb_t lo = p.mantissa << bit;
    const limb_t hi = bit != 0 ? (p.mantissa >> (kLimbBits - bit)) : 0;

    if (!limbs::add_shifted(entry(c, i), c.limbs, lo, hi, off, subtract)) {
        return;
    }

    // Signed overflow. Undo it -- two's complement addition is exactly
    // invertible -- then widen and redo. One extra limb always suffices: with
    // n limbs both operands are bounded by 2^(64n-1), so the sum needs at most
    // 64n + 1 bits. That makes this a single retry, never a loop.
    limbs::add_shifted(entry(c, i), c.limbs, lo, hi, off, !subtract);
    ensure_limbs(c, c.limbs + 1);

    const bool still_overflows =
        limbs::add_shifted(entry(c, i), c.limbs, lo, hi, off, subtract);
    assert(!still_overflows && "one extra limb must absorb the carry");
    (void)still_overflows;
}

void ColumnBlockMatrix::add(std::size_t i, std::size_t j, double v)
{
    accumulate(i, j, v, false, 0);
}

void ColumnBlockMatrix::sub(std::size_t i, std::size_t j, double v)
{
    accumulate(i, j, v, true, 0);
}

void ColumnBlockMatrix::add_matrix(const double* b, std::size_t row_stride)
{
    add_matrix_scaled_pow2(b, 0, row_stride);
}

void ColumnBlockMatrix::add_matrix_scaled_pow2(const double* b, int log2_scale,
                                               std::size_t row_stride)
{
    const std::size_t stride = row_stride ? row_stride : cols_;
    // Column-major traversal: each column's scale is resolved once and the
    // block it owns stays hot in cache.
    for (std::size_t j = 0; j < cols_; ++j) {
        for (std::size_t i = 0; i < rows_; ++i) {
            accumulate(i, j, b[i * stride + j], false, log2_scale);
        }
    }
}

void ColumnBlockMatrix::set_zero()
{
    for (auto& c : cols_state_) std::fill(c.data.begin(), c.data.end(), 0);
}

std::vector<limb_t> ColumnBlockMatrix::magnitude(std::size_t i, std::size_t j,
                                                 bool* negative) const
{
    check_index(i, j);
    const Column& c = cols_state_[j];
    const limb_t* src = entry(c, i);
    const bool neg = limbs::is_negative(src, c.limbs);
    if (negative) *negative = neg;

    // One extra limb keeps the magnitude unambiguously non-negative, so the
    // sign-aware helpers can be reused on it.
    std::vector<limb_t> mag(c.limbs + 1, 0);
    if (neg) {
        limbs::negate_into(src, c.limbs, mag.data());
        mag[c.limbs] = 0;
    } else {
        std::copy(src, src + c.limbs, mag.begin());
    }
    return mag;
}

bool ColumnBlockMatrix::is_zero(std::size_t i, std::size_t j) const
{
    check_index(i, j);
    const Column& c = cols_state_[j];
    return limbs::is_zero(entry(c, i), c.limbs);
}

double ColumnBlockMatrix::to_double(std::size_t i, std::size_t j) const
{
    bool neg = false;
    const std::vector<limb_t> mag = magnitude(i, j, &neg);
    const std::size_t b = limbs::bit_length(mag.data(), mag.size());
    if (b == 0) return 0.0;

    const long long exp = cols_state_[j].exponent;
    const long long top = static_cast<long long>(b) - 1 + exp;  // 2^top <= |x|

    // Normal range: keep 53 significant bits. Subnormal range: round directly
    // to a multiple of 2^-1074, which avoids the double rounding a two-step
    // (round to 53 bits, then let ldexp round again) would introduce.
    const long long drop =
        (top < -1022) ? (-1074 - exp) : (static_cast<long long>(b) - 53);

    double r;
    if (drop <= 0) {
        const std::uint64_t m = limbs::extract_u64(mag.data(), mag.size(), 0);
        r = std::ldexp(static_cast<double>(m), static_cast<int>(exp));
    } else {
        const std::size_t d = static_cast<std::size_t>(drop);
        const bool round_bit = limbs::get_bit(mag.data(), mag.size(), d - 1);
        const bool sticky =
            limbs::any_bits_below(mag.data(), mag.size(), d - 1);
        std::uint64_t m = limbs::extract_u64(mag.data(), mag.size(), d);
        if (round_bit && (sticky || (m & 1) != 0)) ++m;
        r = std::ldexp(static_cast<double>(m), static_cast<int>(exp + drop));
    }
    return neg ? -r : r;
}

void ColumnBlockMatrix::to_matrix(double* out, std::size_t row_stride) const
{
    const std::size_t stride = row_stride ? row_stride : cols_;
    for (std::size_t j = 0; j < cols_; ++j) {
        for (std::size_t i = 0; i < rows_; ++i) {
            out[i * stride + j] = to_double(i, j);
        }
    }
}

bool ColumnBlockMatrix::is_exactly_representable(std::size_t i,
                                                 std::size_t j) const
{
    bool neg = false;
    const std::vector<limb_t> mag = magnitude(i, j, &neg);
    const std::size_t b = limbs::bit_length(mag.data(), mag.size());
    if (b == 0) return true;

    const long long exp = cols_state_[j].exponent;
    const std::size_t tz = trailing_zeros(mag);
    if (b - tz > 53) return false;                                 // too wide
    if (exp + static_cast<long long>(b) - 1 > 1023) return false;  // overflow
    if (exp + static_cast<long long>(tz) < -1074) return false;    // underflow
    return true;
}

std::string ColumnBlockMatrix::to_exact_decimal(std::size_t i,
                                                std::size_t j) const
{
    bool neg = false;
    std::vector<limb_t> mag = magnitude(i, j, &neg);
    if (limbs::bit_length(mag.data(), mag.size()) == 0) return "0";

    const long long exp = cols_state_[j].exponent;
    std::string digits;
    std::size_t frac_digits = 0;

    if (exp >= 0) {
        const std::size_t bits = limbs::bit_length(mag.data(), mag.size()) +
                                 static_cast<std::size_t>(exp);
        std::vector<limb_t> shifted(limbs_for_bits(bits) + 1, 0);
        limbs::shift_left(mag.data(), mag.size(), static_cast<unsigned>(exp),
                          shifted.data(), shifted.size());
        digits = limbs::magnitude_to_decimal(std::move(shifted));
    } else {
        // v * 2^-k == (v * 5^k) / 10^k, so the expansion is finite. Multiply
        // by 5^k a limb-sized power at a time.
        constexpr std::size_t kStride = 27;
        const std::size_t k = static_cast<std::size_t>(-exp);
        frac_digits = k;
        for (std::size_t done = 0; done < k; done += kStride) {
            limbs::mul_small(mag, pow5(std::min(kStride, k - done)));
        }
        digits = limbs::magnitude_to_decimal(std::move(mag));
    }

    std::string out;
    if (frac_digits == 0) {
        out = digits;
    } else {
        if (digits.size() <= frac_digits) {
            digits.insert(digits.begin(), frac_digits + 1 - digits.size(), '0');
        }
        out = digits.substr(0, digits.size() - frac_digits) + "." +
              digits.substr(digits.size() - frac_digits);
        out.erase(out.find_last_not_of('0') + 1);
        if (!out.empty() && out.back() == '.') out.pop_back();
    }
    return neg ? "-" + out : out;
}

void ColumnBlockMatrix::reserve_column(std::size_t j, int exponent,
                                       std::size_t bits)
{
    if (j >= cols_) throw std::out_of_range("cbfp: column index out of range");
    Column& c = cols_state_[j];
    if (!c.initialized) {
        c.exponent = exponent;
        c.initialized = true;
    } else if (exponent < c.exponent) {
        rescale(c, exponent);
    }
    ensure_limbs(c, limbs_for_bits(bits));
}

void ColumnBlockMatrix::reserve_for(const double* b, std::size_t count,
                                    std::size_t row_stride)
{
    const std::size_t stride = row_stride ? row_stride : cols_;
    const std::size_t headroom = ceil_log2(std::max<std::size_t>(count, 1)) + 1;

    for (std::size_t j = 0; j < cols_; ++j) {
        bool any = false;
        long long low = 0, high = 0;
        for (std::size_t i = 0; i < rows_; ++i) {
            const double v = b[i * stride + j];
            if (!std::isfinite(v)) continue;
            const DoubleParts p = decompose(v);
            if (p.mantissa == 0) continue;

            const long long lsb = p.exponent;
            const long long msb =
                p.exponent +
                static_cast<long long>(kLimbBits - __builtin_clzll(p.mantissa));
            if (!any) {
                low = lsb;
                high = msb;
                any = true;
            } else {
                low = std::min(low, lsb);
                high = std::max(high, msb);
            }
        }
        if (!any) continue;
        reserve_column(j, static_cast<int>(low),
                       static_cast<std::size_t>(high - low) + headroom);
    }
}

std::size_t ColumnBlockMatrix::memory_bytes() const
{
    std::size_t total = sizeof(*this);
    for (const auto& c : cols_state_) {
        total += sizeof(Column) + c.data.capacity() * sizeof(limb_t);
    }
    return total;
}

std::string ColumnBlockMatrix::describe() const
{
    std::ostringstream os;
    os << rows_ << " x " << cols_ << " column-block fixed-point matrix ("
       << memory_bytes() << " bytes)\n";
    for (std::size_t j = 0; j < cols_; ++j) {
        const Column& c = cols_state_[j];
        os << "  col " << j << ": exponent 2^" << c.exponent << ", width "
           << c.limbs * kLimbBits << " bits";
        if (!c.initialized) os << " (unset)";
        os << "\n";
    }
    return os.str();
}

}  // namespace cbfp
