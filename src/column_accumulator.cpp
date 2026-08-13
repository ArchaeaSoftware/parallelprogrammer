#include "cbfp/column_accumulator.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "kernels.hpp"

namespace cbfp {

using limbs::kLimbBits;
using limbs::limb_t;

namespace {

std::size_t
limbs_for_bits(std::size_t bits)
{
    return (bits + kLimbBits - 1) / kLimbBits;
}

// Number of trailing zero bits of a non-negative magnitude (0 if zero).
std::size_t
trailing_zeros(const std::vector<limb_t> &v)
{
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (0 != v[i]) return i * kLimbBits + __builtin_ctzll(v[i]);
    }
    return 0;
}

// Significant bits of a 64-bit quantity. Deliberately not expressed in terms
// of kLimbBits: the callers pass a std::size_t and a double's mantissa,
// neither of which is a limb.
std::size_t
bit_width_u64(std::uint64_t v)
{
    return 0 == v ? 0 : 64 - static_cast<std::size_t>(__builtin_clzll(v));
}

std::size_t
ceil_log2(std::size_t n)
{
    if (n <= 1) return 0;
    return bit_width_u64(n - 1);
}

// 5^k for k <= 27, the largest power of five that fits in a limb.
limb_t
pow5(std::size_t k)
{
    assert(k <= 27);
    limb_t p = 1;
    for (std::size_t n = 0; n < k; ++n) p *= 5;
    return p;
}

constexpr long long kExponentLimit = 1 << 24;

}  // namespace

const char *
active_kernel()
{
    return kernels::accumulate_name();
}

DoubleParts
decompose(double v)
{
    // Read the IEEE-754 fields directly. frexp/ldexp would do the same job but
    // are real libm calls; this is branch-free apart from the odd-normalizing
    // shift and measures about 9x faster.
    //
    //   normal    (biased != 0): m = 2^52 | frac, e = biased - 1075
    //   subnormal (biased == 0): m = frac,        e = -1074
    //
    // Both collapse to e = max(biased, 1) - 1075.
    std::uint64_t bits;
    std::memcpy(&bits, &v, sizeof bits);

    const std::uint64_t biased = (bits >> 52) & 0x7FF;
    const std::uint64_t frac = bits & ((std::uint64_t{1} << 52) - 1);

    std::uint64_t m = frac;
    if (0 != biased) m |= std::uint64_t{1} << 52;
    int e = static_cast<int>(std::max<std::uint64_t>(biased, 1)) - 1075;

    // Normalizing the mantissa to odd makes the exponent the value's true ulp,
    // which keeps columns of coarse values (integers, say) narrow.
    if (0 != m) {
        const int tz = __builtin_ctzll(m);
        m >>= tz;
        e += tz;
    }

    DoubleParts p;
    p.mantissa = m;
    p.negative = (bits >> 63) != 0;
    p.exponent = e;
    return p;
}

ColumnBlockMatrix::ColumnBlockMatrix(std::size_t rows, std::size_t cols)
    : rows_(rows), cols_(cols), cols_state_(cols)
{
    for (auto &c : cols_state_) {
        c.limbs.emplace_back(rows_, 0);
        rebuild_bases(c);
    }
    column_buffer_.resize(rows_);
}

void
ColumnBlockMatrix::check_index(std::size_t i, std::size_t j) const
{
    if (i >= rows_ || j >= cols_) {
        throw std::out_of_range("cbfp: matrix index out of range");
    }
}

void
ColumnBlockMatrix::rebuild_bases(Column &c)
{
    c.bases.resize(c.limbs.size());
    for (std::size_t k = 0; k < c.limbs.size(); ++k) {
        c.bases[k] = c.limbs[k].data();
    }
}

void
ColumnBlockMatrix::ensure_limb_count(Column &c, std::size_t needed)
{
    if (c.limbs.size() >= needed) return;

    // Appending never touches the limb arrays already there -- the whole point
    // of one allocation per limb position. Each new one is sign-filled from
    // the limb below it so negative entries stay negative.
    while (c.limbs.size() < needed) {
        const std::size_t index = c.limbs.size();
        c.limbs.emplace_back(rows_, index);
        kernels::sign_fill(c.limbs[index].data(), c.limbs[index - 1].data(),
                           rows_);
    }
    rebuild_bases(c);
}

void
ColumnBlockMatrix::fit_column(Column &c)
{
    const std::size_t bits =
        c.max_addend_bits + ceil_log2(c.add_count + 1) + 1;  // +1 for the sign
    ensure_limb_count(c, std::max<std::size_t>(limbs_for_bits(bits), 1));
}

bool
ColumnBlockMatrix::column_is_zero(const Column &c) const
{
    for (const auto &lc : c.limbs) {
        const limb_t *p = lc.data();
        for (std::size_t i = 0; i < rows_; ++i) {
            if (0 != p[i]) return false;
        }
    }
    return true;
}

void
ColumnBlockMatrix::rescale(Column &c, int new_exponent)
{
    if (new_exponent >= c.exponent) return;
    const unsigned shift = static_cast<unsigned>(
        static_cast<long long>(c.exponent) - new_exponent);

    // An all-zero column carries no information, so its scale moves freely.
    if (column_is_zero(c)) {
        c.exponent = new_exponent;
        c.max_addend_bits = 0;
        c.add_count = 0;
        return;
    }

    const std::size_t widened =
        c.max_addend_bits + shift + ceil_log2(c.add_count + 1) + 1;
    const std::size_t ndst =
        std::max(limbs_for_bits(widened), c.limbs.size() + shift / kLimbBits);

    std::vector<LimbColumn> fresh;
    fresh.reserve(ndst);
    for (std::size_t k = 0; k < ndst; ++k) fresh.emplace_back(rows_, k);

    std::vector<limb_t *> fresh_bases(ndst);
    for (std::size_t k = 0; k < ndst; ++k) fresh_bases[k] = fresh[k].data();

    kernels::shift_left(fresh_bases.data(), ndst, c.bases.data(),
                        c.limbs.size(), rows_, shift);

    c.limbs = std::move(fresh);
    c.bases = std::move(fresh_bases);
    c.max_addend_bits += shift;
    c.exponent = new_exponent;
}

void
ColumnBlockMatrix::accumulate(std::size_t i, std::size_t j, double v,
                              bool negate, int log2_scale)
{
    check_index(i, j);
    if (!std::isfinite(v)) {
        throw std::domain_error("cbfp: cannot accumulate a non-finite value");
    }

    const DoubleParts p = decompose(v);
    if (0 == p.mantissa) return;

    const long long e64 = static_cast<long long>(p.exponent) + log2_scale;
    if (e64 < -kExponentLimit || e64 > kExponentLimit) {
        throw std::domain_error("cbfp: log2_scale puts the value out of range");
    }
    const int e = static_cast<int>(e64);

    Column &c = cols_state_[j];
    if (!c.initialized) {
        c.exponent = e;
        c.initialized = true;
    }
    if (e < c.exponent) rescale(c, e);

    const std::size_t shift = static_cast<std::size_t>(e - c.exponent);
    c.max_addend_bits =
        std::max(c.max_addend_bits, shift + bit_width_u64(p.mantissa));
    ++c.add_count;
    fit_column(c);

    kernels::accumulate_one(c.bases.data(), c.limbs.size(), i, p.mantissa,
                            shift, p.negative != negate);
}

void
ColumnBlockMatrix::add(std::size_t i, std::size_t j, double v)
{
    accumulate(i, j, v, false, 0);
}

void
ColumnBlockMatrix::sub(std::size_t i, std::size_t j, double v)
{
    accumulate(i, j, v, true, 0);
}

void
ColumnBlockMatrix::add_matrix(const double *b, std::size_t row_stride)
{
    add_matrix_scaled_pow2(b, 0, row_stride);
}

void
ColumnBlockMatrix::add_matrix_scaled_pow2(const double *b, int log2_scale,
                                          std::size_t row_stride)
{
    accumulate_columns(b, 1, row_stride ? row_stride : cols_, log2_scale);
}

void
ColumnBlockMatrix::add_matrix_col_major(const double *b, std::size_t col_stride)
{
    add_matrix_col_major_scaled_pow2(b, 0, col_stride);
}

void
ColumnBlockMatrix::add_matrix_col_major_scaled_pow2(const double *b,
                                                    int log2_scale,
                                                    std::size_t col_stride)
{
    accumulate_columns(b, col_stride ? col_stride : rows_, 1, log2_scale);
}

void
ColumnBlockMatrix::accumulate_columns(const double *b, std::size_t column_step,
                                      std::size_t row_step, int log2_scale)
{
    if (0 == rows_ || 0 == cols_) return;

    for (std::size_t j = 0; j < cols_; ++j) {
        // Both passes want the column contiguous. When the caller's layout
        // already provides that, use it in place; otherwise stage it, because
        // a strided vector gather costs more than the work it feeds.
        const double *column = b + j * column_step;
        if (1 != row_step) {
            for (std::size_t i = 0; i < rows_; ++i) {
                column_buffer_[i] = column[i * row_step];
            }
            column = column_buffer_.data();
        }
        accumulate_column(j, column, log2_scale);
    }
}

void
ColumnBlockMatrix::add_column(std::size_t j, const double *v)
{
    add_column_scaled_pow2(j, v, 0);
}

void
ColumnBlockMatrix::add_column_scaled_pow2(std::size_t j, const double *v,
                                          int log2_scale)
{
    if (j >= cols_) throw std::out_of_range("cbfp: column index out of range");
    if (0 == rows_) return;
    accumulate_column(j, v, log2_scale);
}

void
ColumnBlockMatrix::accumulate_column(std::size_t j, const double *column,
                                     int log2_scale)
{
    // First pass learns only the exponent range, because the rescale and widen
    // decisions have to be made before any value can be added.
    // Once the column has a scale, a lower bound on the incoming exponents is
    // enough to rule out a rescale, and the survey can skip the significand.
    Column &c = cols_state_[j];
    // A column with no scale yet takes whatever the survey returns as its
    // exponent, so there the exact value is always required -- a floor nothing
    // can clear.
    const long long floor_exponent =
        c.initialized ? static_cast<long long>(c.exponent) - log2_scale
                      : std::numeric_limits<long long>::max();
    const kernels::Survey sc = kernels::survey()(column, rows_, floor_exponent);
    if (sc.nonfinite) {
        throw std::domain_error("cbfp: cannot accumulate a non-finite value");
    }
    if (!sc.any) return;

    const long long min_exponent = sc.min_exponent + log2_scale;
    const long long max_top = sc.max_top + log2_scale;
    if (min_exponent < -kExponentLimit || max_top > kExponentLimit) {
        throw std::domain_error("cbfp: log2_scale puts the value out of range");
    }

    if (!c.initialized) {
        c.exponent = static_cast<int>(min_exponent);
        c.initialized = true;
    }
    if (min_exponent < c.exponent) rescale(c, static_cast<int>(min_exponent));

    c.max_addend_bits = std::max(
        c.max_addend_bits, static_cast<std::size_t>(max_top - c.exponent));
    ++c.add_count;
    fit_column(c);

    // Second pass fuses decomposition into the add, so no decomposed form is
    // ever written to memory. Folding the scale into the column exponent keeps
    // the kernel free of it.
    kernels::accumulate()(
        c.bases.data(), c.limbs.size(), column, rows_,
        static_cast<std::int32_t>(c.exponent - log2_scale),
        static_cast<std::size_t>(min_exponent - c.exponent) / kLimbBits);
}

void
ColumnBlockMatrix::set_zero()
{
    for (auto &c : cols_state_) {
        for (auto &lc : c.limbs) {
            std::memset(lc.data(), 0, padded_rows(rows_) * sizeof(limb_t));
        }
        c.max_addend_bits = 0;
        c.add_count = 0;
    }
}

std::vector<limb_t>
ColumnBlockMatrix::magnitude(std::size_t i, std::size_t j, bool *negative) const
{
    check_index(i, j);
    const Column &c = cols_state_[j];
    const std::size_t n = c.limbs.size();

    std::vector<limb_t> value(n);
    for (std::size_t k = 0; k < n; ++k) value[k] = c.limbs[k].data()[i];

    const bool neg = limbs::is_negative(value.data(), n);
    if (negative) *negative = neg;

    // One extra limb keeps the magnitude unambiguously non-negative, so the
    // sign-aware helpers can be reused on it.
    std::vector<limb_t> mag(n + 1, 0);
    if (neg) {
        limbs::negate_into(value.data(), n, mag.data());
    } else {
        std::copy(value.begin(), value.end(), mag.begin());
    }
    return mag;
}

std::vector<limb_t>
ColumnBlockMatrix::entry_limbs(std::size_t i, std::size_t j) const
{
    check_index(i, j);
    const Column &c = cols_state_[j];
    std::vector<limb_t> v(c.limbs.size());
    for (std::size_t k = 0; k < c.limbs.size(); ++k) {
        v[k] = c.limbs[k].data()[i];
    }
    return v;
}

bool
ColumnBlockMatrix::is_zero(std::size_t i, std::size_t j) const
{
    check_index(i, j);
    const Column &c = cols_state_[j];
    for (const auto &lc : c.limbs) {
        if (0 != lc.data()[i]) return false;
    }
    return true;
}

double
ColumnBlockMatrix::to_double(std::size_t i, std::size_t j) const
{
    bool neg = false;
    const std::vector<limb_t> mag = magnitude(i, j, &neg);
    const std::size_t b = limbs::bit_length(mag.data(), mag.size());
    if (0 == b) return 0.0;

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
        if (round_bit && (sticky || 0 != (m & 1))) ++m;
        r = std::ldexp(static_cast<double>(m), static_cast<int>(exp + drop));
    }
    return neg ? -r : r;
}

void
ColumnBlockMatrix::to_matrix(double *out, std::size_t row_stride) const
{
    const std::size_t stride = row_stride ? row_stride : cols_;
    for (std::size_t j = 0; j < cols_; ++j) {
        for (std::size_t i = 0; i < rows_; ++i) {
            out[i * stride + j] = to_double(i, j);
        }
    }
}

bool
ColumnBlockMatrix::is_exactly_representable(std::size_t i, std::size_t j) const
{
    bool neg = false;
    const std::vector<limb_t> mag = magnitude(i, j, &neg);
    const std::size_t b = limbs::bit_length(mag.data(), mag.size());
    if (0 == b) return true;

    const long long exp = cols_state_[j].exponent;
    const std::size_t tz = trailing_zeros(mag);
    if (b - tz > 53) return false;                                 // too wide
    if (exp + static_cast<long long>(b) - 1 > 1023) return false;  // overflow
    if (exp + static_cast<long long>(tz) < -1074) return false;    // underflow
    return true;
}

std::string
ColumnBlockMatrix::to_exact_decimal(std::size_t i, std::size_t j) const
{
    bool neg = false;
    std::vector<limb_t> mag = magnitude(i, j, &neg);
    if (0 == limbs::bit_length(mag.data(), mag.size())) return "0";

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
    if (0 == frac_digits) {
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

void
ColumnBlockMatrix::reserve_column(std::size_t j, int exponent, std::size_t bits)
{
    if (j >= cols_) throw std::out_of_range("cbfp: column index out of range");
    Column &c = cols_state_[j];
    if (!c.initialized) {
        c.exponent = exponent;
        c.initialized = true;
    } else if (exponent < c.exponent) {
        rescale(c, exponent);
    }
    ensure_limb_count(c, std::max<std::size_t>(limbs_for_bits(bits), 1));
}

void
ColumnBlockMatrix::reserve_for(const double *b, std::size_t count,
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
            if (0 == p.mantissa) continue;

            const long long lsb = p.exponent;
            const long long msb =
                p.exponent + static_cast<long long>(bit_width_u64(p.mantissa));
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

std::size_t
ColumnBlockMatrix::memory_bytes() const
{
    std::size_t total = sizeof(*this);
    const std::size_t per_limb = padded_rows(rows_) * sizeof(limb_t) + 512;
    for (const auto &c : cols_state_) {
        total += sizeof(Column) + c.limbs.size() * per_limb +
                 c.bases.capacity() * sizeof(limb_t *);
    }
    return total;
}

std::string
ColumnBlockMatrix::describe() const
{
    std::ostringstream os;
    os << rows_ << " x " << cols_ << " column-block fixed-point matrix ("
       << memory_bytes() << " bytes, " << active_kernel() << " kernel)\n";
    for (std::size_t j = 0; j < cols_; ++j) {
        const Column &c = cols_state_[j];
        os << "  col " << j << ": exponent 2^" << c.exponent << ", width "
           << c.limbs.size() * kLimbBits << " bits";
        if (!c.initialized) os << " (unset)";
        os << "\n";
    }
    return os.str();
}

}  // namespace cbfp
