#include "truesum/accumulation_matrix.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <condition_variable>
#include <cstring>
#include <functional>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>

#include "kernels.hpp"
#include "thread_pool.hpp"

namespace truesum {

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
    //   denormal (biased == 0): m = frac,        e = -1074
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

AccumulationMatrix::AccumulationMatrix(std::size_t rows, std::size_t cols)
    : rows_(rows), cols_(cols), cols_state_(cols)
{
    for (std::size_t j = 0; j < cols_; ++j) {
        Column &c = cols_state_[j];
        c.rows = rows_;
        c.first_row = 0;
        c.limbs.emplace_back(c.rows, 0);
        rebuild_bases(c);
    }
    column_buffers_.resize(1);
    column_buffers_[0].resize(rows_);
}

AccumulationMatrix::AccumulationMatrix(std::size_t n, Uplo uplo)
    : rows_(n), cols_(n), symmetric_(true), uplo_(uplo), cols_state_(n)
{
    for (std::size_t j = 0; j < cols_; ++j) {
        Column &c = cols_state_[j];
        c.rows = Uplo::Lower == uplo_ ? n - j : j + 1;
        c.first_row = Uplo::Lower == uplo_ ? j : 0;
        c.limbs.emplace_back(c.rows, 0);
        rebuild_bases(c);
    }
    // Sized for the longest column, so one buffer serves every column.
    column_buffers_.resize(1);
    column_buffers_[0].resize(rows_);
}

std::size_t
AccumulationMatrix::column_rows(std::size_t j) const
{
    if (j >= cols_)
        throw std::out_of_range("truesum: column index out of range");
    return cols_state_[j].rows;
}

std::size_t
AccumulationMatrix::column_first_row(std::size_t j) const
{
    if (j >= cols_)
        throw std::out_of_range("truesum: column index out of range");
    return cols_state_[j].first_row;
}

// Lower stores i >= j, so (i, j) and (j, i) both land in column min(i, j) at
// slot |i - j|. Upper stores i <= j, so both land in column max(i, j) at slot
// min(i, j). Neither needs a branch on which side of the diagonal it started.
void
AccumulationMatrix::locate(std::size_t &col, std::size_t &slot, std::size_t i,
                           std::size_t j) const
{
    if (!symmetric_) {
        col = j;
        slot = i;
        return;
    }
    if (Uplo::Lower == uplo_) {
        col = std::min(i, j);
        slot = i < j ? j - i : i - j;
    } else {
        col = std::max(i, j);
        slot = std::min(i, j);
    }
}

AccumulationMatrix::~AccumulationMatrix() = default;
AccumulationMatrix::AccumulationMatrix(AccumulationMatrix &&) noexcept =
    default;
AccumulationMatrix &
AccumulationMatrix::operator=(AccumulationMatrix &&) noexcept = default;

void
AccumulationMatrix::set_threads(unsigned n)
{
    if (n < 1) n = 1;
    if (n == threads()) return;
    pool_.reset();
    column_buffers_.assign(n, std::vector<double>(rows_));
    if (n > 1) pool_.reset(new detail::ThreadPool(n));
}

unsigned
AccumulationMatrix::threads() const
{
    return static_cast<unsigned>(column_buffers_.size());
}

AccumulationMatrix::AccumulationMatrix(
    std::size_t rows, std::size_t cols,
    const std::vector<std::vector<Survey>> &surveys)
    : AccumulationMatrix(rows, cols)
{
    reserve_from_surveys(surveys);
}

AccumulationMatrix::AccumulationMatrix(
    std::size_t n, Uplo uplo, const std::vector<std::vector<Survey>> &surveys)
    : AccumulationMatrix(n, uplo)
{
    reserve_from_surveys(surveys);
}

void
AccumulationMatrix::reserve_from_surveys(
    const std::vector<std::vector<Survey>> &surveys)
{
    if (surveys.empty()) {
        throw std::invalid_argument(
            "truesum: pre-sizing needs the surveys of at least one matrix");
    }
    for (const auto &one : surveys) {
        if (one.size() != cols_) {
            throw std::invalid_argument(
                "truesum: each survey must have one entry per column");
        }
    }

    // One reservation per column: the minimum of the minima for the exponent,
    // the maximum of the maxima for the top, and the count from the outer
    // size. Reducing to an aggregate is what makes submission order
    // irrelevant -- any matrix inside these extents fits the reservation, so
    // the accumulation matrix never needs to know which one it is being handed.
    const std::size_t headroom = ceil_log2(surveys.size() + 1) + 1;
    for (std::size_t j = 0; j < cols_; ++j) {
        bool any = false;
        int low = 0, high = 0;
        for (const auto &one : surveys) {
            if (!one[j].any) continue;
            if (!any) {
                low = one[j].min_exponent;
                high = one[j].max_top;
                any = true;
            } else {
                low = std::min(low, one[j].min_exponent);
                high = std::max(high, one[j].max_top);
            }
        }
        if (!any) continue;
        reserve_column(j, low,
                       static_cast<std::size_t>(high - low) + headroom);
    }

    presized_ = true;
    declared_matrices_ = surveys.size();
}

// A batch reached outside what its column was sized for. With surveys supplied
// by a producer that means the metadata was wrong; with surveys computed here
// it means the survey and the accumulate disagree, which is a bug. Either way
// the sums are already wrong, so this reports rather than recovers.
void
AccumulationMatrix::report_contradiction(std::size_t j, unsigned flags) const
{
    std::ostringstream os;
    os << "truesum: column " << j << " received values it was not sized for:";
    if (0 != (flags & kernels::kBadNonFinite)) os << " a non-finite value;";
    if (0 != (flags & kernels::kBadExponent)) {
        os << " an exponent below the column's;";
    }
    if (0 != (flags & kernels::kBadWidth)) os << " an addend past its width;";
    if (0 != (flags & kernels::kBadOverflow)) os << " a sum past its width;";
    os << " the accumulation matrix is no longer consistent";
    throw std::runtime_error(os.str());
}

void
AccumulationMatrix::check_index(std::size_t i, std::size_t j) const
{
    if (i >= rows_ || j >= cols_) {
        throw std::out_of_range("truesum: matrix index out of range");
    }
}

void
AccumulationMatrix::rebuild_bases(Column &c)
{
    c.bases.resize(c.limbs.size());
    for (std::size_t k = 0; k < c.limbs.size(); ++k) {
        c.bases[k] = c.limbs[k].data();
    }
}

void
AccumulationMatrix::ensure_limb_count(Column &c, std::size_t needed)
{
    if (c.limbs.size() >= needed) return;

    // Appending never touches the limb arrays already there -- the whole point
    // of one allocation per limb position. Each new one is sign-filled from
    // the limb below it so negative entries stay negative.
    while (c.limbs.size() < needed) {
        const std::size_t index = c.limbs.size();
        c.limbs.emplace_back(c.rows, index);
        kernels::sign_fill(c.limbs[index].data(), c.limbs[index - 1].data(),
                           c.rows);
    }
    rebuild_bases(c);
}

void
AccumulationMatrix::fit_column(Column &c)
{
    const std::size_t bits =
        c.max_addend_bits + ceil_log2(c.add_count + 1) + 1;  // +1 for the sign
    ensure_limb_count(c, std::max<std::size_t>(limbs_for_bits(bits), 1));
}

bool
AccumulationMatrix::column_is_zero(const Column &c) const
{
    for (const auto &lc : c.limbs) {
        const limb_t *p = lc.data();
        for (std::size_t i = 0; i < c.rows; ++i) {
            if (0 != p[i]) return false;
        }
    }
    return true;
}

void
AccumulationMatrix::rescale(Column &c, int new_exponent)
{
    if (new_exponent >= c.exponent) return;
    const unsigned shift = static_cast<unsigned>(
        static_cast<long long>(c.exponent) - new_exponent);

    // An all-zero column holds no information, so its scale moves freely.
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
    for (std::size_t k = 0; k < ndst; ++k) fresh.emplace_back(c.rows, k);

    std::vector<limb_t *> fresh_bases(ndst);
    for (std::size_t k = 0; k < ndst; ++k) fresh_bases[k] = fresh[k].data();

    kernels::shift_left(fresh_bases.data(), ndst, c.bases.data(),
                        c.limbs.size(), c.rows, shift);

    c.limbs = std::move(fresh);
    c.bases = std::move(fresh_bases);
    c.max_addend_bits += shift;
    c.exponent = new_exponent;
}

void
AccumulationMatrix::add_matrix(const double *b, std::size_t row_stride)
{
    add_matrix_scaled_pow2(b, 0, row_stride);
}

void
AccumulationMatrix::add_matrix_scaled_pow2(const double *b, int log2_scale,
                                           std::size_t row_stride)
{
    accumulate_columns(b, 1, row_stride ? row_stride : cols_, log2_scale);
}

void
AccumulationMatrix::add_matrix_col_major(const double *b,
                                         std::size_t col_stride)
{
    add_matrix_col_major_scaled_pow2(b, 0, col_stride);
}

void
AccumulationMatrix::add_matrix_col_major_scaled_pow2(const double *b,
                                                     int log2_scale,
                                                     std::size_t col_stride)
{
    accumulate_columns(b, col_stride ? col_stride : rows_, 1, log2_scale);
}

void
AccumulationMatrix::accumulate_column_range(const double *b,
                                            std::size_t column_step,
                                            std::size_t row_step,
                                            int log2_scale, std::size_t begin,
                                            std::size_t end, unsigned slot)
{
    std::vector<double> &staging = column_buffers_[slot];
    for (std::size_t j = begin; j < end; ++j) {
        // Both passes want the column contiguous. When the caller's layout
        // already provides that, use it in place; otherwise stage it, because
        // a strided vector gather costs more than the work it feeds.
        const Column &c = cols_state_[j];
        const double *column = b + j * column_step + c.first_row * row_step;
        if (1 != row_step) {
            for (std::size_t i = 0; i < c.rows; ++i) {
                staging[i] = column[i * row_step];
            }
            column = staging.data();
        }
        accumulate_column(j, column, log2_scale);
    }
}

// A triangular column's work runs from n entries down to 1, so handing every
// worker the same number of columns hands one of them most of the matrix.
// Assign column j to the slot its cumulative entry count falls in: the count
// rises monotonically, so each slot still gets a contiguous block and a
// worker's columns stay near one another in memory.
void
AccumulationMatrix::partition_columns(std::size_t &begin, std::size_t &end,
                                      unsigned slot) const
{
    const unsigned n = threads();
    if (!symmetric_) {
        pool_->partition(begin, end, cols_, slot);
        return;
    }

    std::size_t total = 0;
    for (const auto &c : cols_state_) total += c.rows;
    if (0 == total) {
        begin = end = 0;
        return;
    }

    begin = end = cols_;
    std::size_t seen = 0;  // entries in the columns before j
    for (std::size_t j = 0; j < cols_; ++j) {
        const unsigned owner = static_cast<unsigned>(seen * n / total);
        if (owner == slot && cols_ == begin) begin = j;
        if (owner > slot) {
            end = j;
            break;
        }
        seen += cols_state_[j].rows;
    }
    if (cols_ == begin) begin = end = cols_;  // this slot drew nothing
}

void
AccumulationMatrix::accumulate_columns(const double *b, std::size_t column_step,
                                       std::size_t row_step, int log2_scale)
{
    if (0 == rows_ || 0 == cols_) return;
    if (presized_ && ++submitted_matrices_ > declared_matrices_) {
        throw std::runtime_error(
            "truesum: more matrices accumulated than were described to the "
            "constructor; the width bound holds for that many and no more");
    }

    const unsigned n = threads();
    if (1 == n || cols_ < 2) {
        accumulate_column_range(b, column_step, row_step, log2_scale, 0, cols_,
                                0);
        return;
    }
    // A contiguous block each, so a worker's columns stay near one another in
    // memory. Columns are independent in storage, so nothing needs locking.
    pool_->run([&](unsigned slot) {
        std::size_t begin = 0, end = 0;
        partition_columns(begin, end, slot);
        accumulate_column_range(b, column_step, row_step, log2_scale, begin,
                                end, slot);
    });
}

void
AccumulationMatrix::add_matrices_col_major(const double *const *b,
                                           std::size_t count,
                                           std::size_t col_stride)
{
    if (0 == count) return;
    if (0 == rows_ || 0 == cols_) return;
    if (presized_) {
        submitted_matrices_ += count;
        if (submitted_matrices_ > declared_matrices_) {
            throw std::runtime_error(
                "truesum: more matrices accumulated than were described to the "
                "constructor; the width bound holds for that many and no more");
        }
    }
    const std::size_t stride = col_stride ? col_stride : rows_;

    const unsigned n = threads();
    if (1 == n || cols_ < 2) {
        fold_column_range(b, count, stride, 0, cols_);
        return;
    }
    pool_->run([&](unsigned slot) {
        std::size_t begin = 0, end = 0;
        partition_columns(begin, end, slot);
        fold_column_range(b, count, stride, begin, end);
    });
}

void
AccumulationMatrix::fold_column_range(const double *const *b, std::size_t count,
                                      std::size_t col_stride, std::size_t begin,
                                      std::size_t end)
{
    // One pointer per matrix, rebuilt per column. Small and on the stack of
    // whichever worker is running, so the columns share nothing.
    std::vector<const double *> columns(count);
    for (std::size_t j = begin; j < end; ++j) {
        const Column &c = cols_state_[j];
        for (std::size_t k = 0; k < count; ++k) {
            columns[k] = b[k] + j * col_stride + c.first_row;
        }
        fold_column(j, columns.data(), count);
    }
}

void
AccumulationMatrix::fold_column(std::size_t j, const double *const *columns,
                                std::size_t count)
{
    Column &c = cols_state_[j];

    // The exponent and width have to be fixed before the fold starts: a
    // rescale partway through would have to re-shift limbs that earlier
    // matrices in this same fold had already been added to. Pre-sizing fixes
    // them in the constructor; without it, surveying all `count` columns first
    // and fitting once fixes them here, which is the same guarantee arrived
    // at later.
    if (!presized_) {
        bool any = false;
        long long low = 0, high = 0;
        for (std::size_t k = 0; k < count; ++k) {
            const kernels::Survey sv = kernels::survey()(
                columns[k], c.rows, std::numeric_limits<long long>::max());
            if (sv.nonfinite) {
                throw std::domain_error(
                    "truesum: cannot accumulate a non-finite value");
            }
            if (!sv.any) continue;
            if (!any) {
                low = sv.min_exponent;
                high = sv.max_top;
                any = true;
            } else {
                low = std::min<long long>(low, sv.min_exponent);
                high = std::max<long long>(high, sv.max_top);
            }
        }
        if (!any) return;

        if (!c.initialized) {
            c.exponent = static_cast<int>(low);
            c.initialized = true;
        }
        if (low < c.exponent) rescale(c, static_cast<int>(low));
        c.max_addend_bits = std::max(
            c.max_addend_bits, static_cast<std::size_t>(high - c.exponent));
        c.add_count += count;
        fit_column(c);
    }

    unsigned flags = 0;
    if (c.limbs.size() <= kernels::kMaxFoldLimbs) {
        kernels::accumulate_fold()(c.bases.data(), c.limbs.size(), columns,
                                   count, c.rows,
                                   static_cast<std::int32_t>(c.exponent),
                                   &flags);
    } else {
        // Too wide to hold a row in registers. One batch at a time is not a
        // fallback so much as the better shape here, since that kernel stops
        // at the first dead carry instead of writing every limb back.
        for (std::size_t k = 0; k < count; ++k) {
            kernels::accumulate()(c.bases.data(), c.limbs.size(), columns[k],
                                  c.rows,
                                  static_cast<std::int32_t>(c.exponent), 0,
                                  &flags);
        }
    }
    if (0 != flags) report_contradiction(j, flags);
}

void
AccumulationMatrix::add_column(std::size_t j, const double *v)
{
    add_column_scaled_pow2(j, v, 0);
}

void
AccumulationMatrix::add_column_scaled_pow2(std::size_t j, const double *v,
                                           int log2_scale)
{
    if (j >= cols_)
        throw std::out_of_range("truesum: column index out of range");
    if (0 == cols_state_[j].rows) return;
    accumulate_column(j, v, log2_scale);
}

void
AccumulationMatrix::accumulate_column(std::size_t j, const double *column,
                                      int log2_scale)
{
    // First pass learns only the exponent range, because the rescale and widen
    // decisions have to be made before any value can be added.
    // Once the column has a scale, a lower bound on the incoming exponents is
    // enough to rule out a rescale, and the survey can skip the significand.
    Column &c = cols_state_[j];
    if (presized_) {
        // Nothing to learn: the extents were supplied, so no rescale or widen
        // can be due and the survey would only confirm what is already known.
        // first_limb is 0 because the column exponent is the aggregate
        // minimum, so every addend starts at or above limb 0.
        unsigned flags = 0;
        kernels::accumulate()(
            c.bases.data(), c.limbs.size(), column, c.rows,
            static_cast<std::int32_t>(c.exponent - log2_scale), 0, &flags);
        if (0 != flags) report_contradiction(j, flags);
        return;
    }
    // A column with no scale yet takes whatever the survey returns as its
    // exponent, so there the exact value is always required -- a floor nothing
    // can clear.
    const long long floor_exponent =
        c.initialized ? static_cast<long long>(c.exponent) - log2_scale
                      : std::numeric_limits<long long>::max();
    const kernels::Survey sc =
        kernels::survey()(column, c.rows, floor_exponent);
    if (sc.nonfinite) {
        throw std::domain_error(
            "truesum: cannot accumulate a non-finite value");
    }
    if (!sc.any) return;

    // Widened before the add, not after: both operands are int, and log2_scale
    // comes from the caller, so summing in int would overflow on the way to the
    // very check below that exists to reject it.
    const long long min_exponent =
        static_cast<long long>(sc.min_exponent) + log2_scale;
    const long long max_top = static_cast<long long>(sc.max_top) + log2_scale;
    if (min_exponent < -kExponentLimit || max_top > kExponentLimit) {
        throw std::domain_error(
            "truesum: log2_scale puts the value out of range");
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
    unsigned flags = 0;
    kernels::accumulate()(
        c.bases.data(), c.limbs.size(), column, c.rows,
        static_cast<std::int32_t>(c.exponent - log2_scale),
        static_cast<std::size_t>(min_exponent - c.exponent) / kLimbBits,
        &flags);
    // Sized from its own survey, this column cannot contradict itself; the
    // check is free and asserts the survey and the accumulate agree.
    if (0 != flags) report_contradiction(j, flags);
}

void
AccumulationMatrix::set_zero()
{
    for (auto &c : cols_state_) {
        for (auto &lc : c.limbs) {
            std::memset(lc.data(), 0, padded_rows(c.rows) * sizeof(limb_t));
        }
        c.max_addend_bits = 0;
        c.add_count = 0;
    }
}

std::vector<limb_t>
AccumulationMatrix::magnitude(bool *negative, std::size_t i,
                              std::size_t j) const
{
    check_index(i, j);
    std::size_t col = 0, slot = 0;
    locate(col, slot, i, j);
    const Column &c = cols_state_[col];
    const std::size_t n = c.limbs.size();

    std::vector<limb_t> value(n);
    for (std::size_t k = 0; k < n; ++k) value[k] = c.limbs[k].data()[slot];

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
AccumulationMatrix::entry_limbs(std::size_t i, std::size_t j) const
{
    check_index(i, j);
    std::size_t col = 0, slot = 0;
    locate(col, slot, i, j);
    const Column &c = cols_state_[col];
    std::vector<limb_t> v(c.limbs.size());
    for (std::size_t k = 0; k < c.limbs.size(); ++k) {
        v[k] = c.limbs[k].data()[slot];
    }
    return v;
}

bool
AccumulationMatrix::is_zero(std::size_t i, std::size_t j) const
{
    check_index(i, j);
    std::size_t col = 0, slot = 0;
    locate(col, slot, i, j);
    const Column &c = cols_state_[col];
    for (const auto &lc : c.limbs) {
        if (0 != lc.data()[slot]) return false;
    }
    return true;
}

namespace {

// A rounded double together with the exact reconstruction of it:
// `value == m * 2^scale`, with `scale >= exp` always. Returning m and scale out
// of the rounding is what lets the residual be an exact limb subtraction --
// `m << (scale - exp)` lands back in the same fixed point the entry is stored
// in, so no second rounding creeps in on the way.
struct Rounded {
    double value;
    std::uint64_t m;
    long long scale;
};

// The correctly rounded double nearest a non-negative magnitude whose bit 0
// has weight 2^exp.
Rounded
round_magnitude(const limb_t *mag, std::size_t n, long long exp)
{
    const std::size_t b = limbs::bit_length(mag, n);
    if (0 == b) return Rounded{0.0, 0, exp};

    const long long top = static_cast<long long>(b) - 1 + exp;  // 2^top <= |x|

    // Normal range: keep 53 significant bits. Denormal range: round directly
    // to a multiple of 2^-1074, which avoids the double rounding a two-step
    // (round to 53 bits, then let ldexp round again) would introduce.
    const long long drop =
        (top < -1022) ? (-1074 - exp) : (static_cast<long long>(b) - 53);

    if (drop <= 0) {
        // Every significant bit fits in the significand, so this is exact and
        // the residual derived from it is zero.
        const std::uint64_t m = limbs::extract_u64(mag, n, 0);
        return Rounded{
            std::ldexp(static_cast<double>(m), static_cast<int>(exp)), m, exp};
    }
    const std::size_t d = static_cast<std::size_t>(drop);
    const bool round_bit = limbs::get_bit(mag, n, d - 1);
    const bool sticky = limbs::any_bits_below(mag, n, d - 1);
    std::uint64_t m = limbs::extract_u64(mag, n, d);
    if (round_bit && (sticky || 0 != (m & 1))) ++m;
    return Rounded{
        std::ldexp(static_cast<double>(m), static_cast<int>(exp + drop)), m,
        exp + drop};
}

// (magnitude - r) at weight 2^exp, rounded to a double. Signed: negative when
// the rounding went up, which is why the caller applies the entry's own sign
// afterwards rather than folding it in here.
double
residual_of(const std::vector<limb_t> &mag, long long exp, const Rounded &r)
{
    // Nothing was discarded: the value is zero, or every bit fit.
    if (0 == r.m || r.scale == exp) return 0.0;
    // Beyond double's range the difference is not representable either.
    if (!std::isfinite(r.value)) return 0.0;

    // Place r.m back at its own weight and subtract, the same 128-bit-at-an-
    // offset decomposition the accumulate kernel uses. `mag` has a spare
    // top limb from magnitude(), so the difference has room for its sign.
    const std::size_t shift = static_cast<std::size_t>(r.scale - exp);
    const std::size_t off = shift / limbs::kLimbBits;
    const unsigned bit = static_cast<unsigned>(shift % limbs::kLimbBits);
    const limb_t lo = 0 == bit ? r.m : r.m << bit;
    const limb_t hi = 0 == bit ? 0 : r.m >> (limbs::kLimbBits - bit);

    std::vector<limb_t> diff = mag;
    limbs::add_shifted(diff.data(), diff.size(), lo, hi, off, true);

    const bool neg = limbs::is_negative(diff.data(), diff.size());
    std::vector<limb_t> d(diff.size());
    if (neg) {
        limbs::negate_into(diff.data(), diff.size(), d.data());
    } else {
        d = diff;
    }
    const double v = round_magnitude(d.data(), d.size(), exp).value;
    return neg ? -v : v;
}

// The correctly rounded double nearest (magnitude * 2^exp) / divisor, for a
// non-negative magnitude and a divisor of at least 1.
//
// A quotient is not a fixed-point value, so it cannot be handed to
// round_magnitude where it lies. Instead the magnitude is shifted left by s
// and divided, which yields floor(magnitude * 2^s / divisor) exactly, and the
// remainder collapses into one sticky bit appended below the quotient. The
// result holds every bit the rounding can look at: the quotient's own bits
// verbatim, and "something nonzero lies below" as the sticky bit. s is chosen
// so the quotient has at least 55 bits, which puts the sticky bit below the
// round bit on both the normal path (53 kept, so at least 3 dropped) and the
// denormal path, where a 56-bit value whose top is under 2^-1022 has its
// lowest bit under 2^-1077.
Rounded
round_quotient(const limb_t *mag, std::size_t n, long long exp,
               std::uint64_t divisor)
{
    const std::size_t bm = limbs::bit_length(mag, n);
    if (0 == bm) return Rounded{0.0, 0, exp};

    // floor(mag * 2^s / divisor) has at least bm + s - bd significant bits.
    const std::size_t bd = bit_width_u64(divisor);
    const std::size_t s = bd + 55 > bm ? bd + 55 - bm : 0;

    std::vector<limb_t> q(limbs_for_bits(bm + s) + 1, 0);
    limbs::shift_left(mag, n, static_cast<unsigned>(s), q.data(), q.size());
    const limb_t rem = limbs::div_small(q, divisor);

    // (q << 1) | sticky, at weight 2^(exp - s - 1).
    std::vector<limb_t> q2(q.size() + 1, 0);
    limbs::shift_left(q.data(), q.size(), 1, q2.data(), q2.size());
    if (0 != rem) q2[0] |= 1;
    return round_magnitude(q2.data(), q2.size(),
                           exp - static_cast<long long>(s) - 1);
}

// (magnitude * 2^exp) / divisor - r, rounded to a double. Signed like
// residual_of, and for the same reason.
//
// The difference is formed as an integer over the same divisor. With t the
// lower of the two weights involved,
//
//     D = magnitude * 2^(exp - t) - r.m * divisor * 2^(r.scale - t)
//
// is an integer, the residual is D * 2^t / divisor, and round_quotient turns
// that into a double with a single rounding.
double
residual_of_quotient(const std::vector<limb_t> &mag, long long exp,
                     std::uint64_t divisor, const Rounded &r)
{
    if (0 == r.m) return 0.0;
    if (!std::isfinite(r.value)) return 0.0;

    const long long t = std::min(exp, r.scale);
    const std::size_t sa = static_cast<std::size_t>(exp - t);
    const std::size_t sb = static_cast<std::size_t>(r.scale - t);

    // r.m * divisor is under 118 bits; the extra limb is room for the sign.
    const std::size_t bm = limbs::bit_length(mag.data(), mag.size());
    const std::size_t n = limbs_for_bits(std::max(bm + sa, sb + 118)) + 1;

    std::vector<limb_t> a(n, 0);
    limbs::shift_left(mag.data(), mag.size(), static_cast<unsigned>(sa),
                      a.data(), n);

    std::vector<limb_t> b(1, r.m);
    limbs::mul_small(b, divisor);
    b.push_back(0);  // keeps it non-negative for the sign-aware shift
    std::vector<limb_t> bs(n, 0);
    limbs::shift_left(b.data(), b.size(), static_cast<unsigned>(sb), bs.data(),
                      n);

    limbs::sub(a.data(), n, bs.data());
    const bool neg = limbs::is_negative(a.data(), n);
    std::vector<limb_t> d(n);
    if (neg) {
        limbs::negate_into(a.data(), n, d.data());
    } else {
        d = a;
    }
    // A negative difference that rounds to zero must not surface as -0.0:
    // the caller applies the entry's sign only to a nonzero residual, and
    // this is the same rule one level down.
    const double v = round_quotient(d.data(), n, t, divisor).value;
    return (neg && 0.0 != v) ? -v : v;
}

}  // namespace

double
AccumulationMatrix::to_double(std::size_t i, std::size_t j) const
{
    return to_double(nullptr, i, j);
}

double
AccumulationMatrix::to_double(double *residual, std::size_t i,
                              std::size_t j) const
{
    bool neg = false;
    const std::vector<limb_t> mag = magnitude(&neg, i, j);
    std::size_t ecol = 0, eslot = 0;
    locate(ecol, eslot, i, j);
    const long long exp = cols_state_[ecol].exponent;
    const Rounded r = round_magnitude(mag.data(), mag.size(), exp);

    if (nullptr != residual) {
        const double lo = residual_of(mag, exp, r);
        // A zero residual is +0.0 whatever the entry's sign. Negating it would
        // hand back -0.0, which reads as a rounding direction that was never
        // taken and compares unequal to the +0.0 an exact readback produces.
        *residual = (neg && 0.0 != lo) ? -lo : lo;
    }
    return neg ? -r.value : r.value;
}

void
AccumulationMatrix::to_matrix(double *out, std::size_t row_stride) const
{
    const std::size_t stride = row_stride ? row_stride : cols_;
    for (std::size_t j = 0; j < cols_; ++j) {
        for (std::size_t i = 0; i < rows_; ++i) {
            out[i * stride + j] = to_double(i, j);
        }
    }
}

void
AccumulationMatrix::to_matrix_with_residual(double *out, double *residual,
                                            std::size_t row_stride) const
{
    const std::size_t stride = row_stride ? row_stride : cols_;
    for (std::size_t j = 0; j < cols_; ++j) {
        for (std::size_t i = 0; i < rows_; ++i) {
            out[i * stride + j] =
                to_double(&residual[i * stride + j], i, j);
        }
    }
}

double
AccumulationMatrix::to_double_mean(std::size_t i, std::size_t j,
                                  std::uint64_t count) const
{
    return to_double_mean(nullptr, i, j, count);
}

double
AccumulationMatrix::to_double_mean(double *residual, std::size_t i,
                                  std::size_t j, std::uint64_t count) const
{
    if (0 == count) throw std::domain_error("truesum: mean over zero values");
    bool neg = false;
    const std::vector<limb_t> mag = magnitude(&neg, i, j);
    std::size_t ecol = 0, eslot = 0;
    locate(ecol, eslot, i, j);
    const long long exp = cols_state_[ecol].exponent;
    const Rounded r = round_quotient(mag.data(), mag.size(), exp, count);

    if (nullptr != residual) {
        const double lo = residual_of_quotient(mag, exp, count, r);
        *residual = (neg && 0.0 != lo) ? -lo : lo;
    }
    return neg ? -r.value : r.value;
}

void
AccumulationMatrix::to_matrix_mean(double *out, std::uint64_t count,
                                  std::size_t row_stride) const
{
    const std::size_t stride = row_stride ? row_stride : cols_;
    for (std::size_t j = 0; j < cols_; ++j) {
        for (std::size_t i = 0; i < rows_; ++i) {
            out[i * stride + j] = to_double_mean(i, j, count);
        }
    }
}

void
AccumulationMatrix::to_matrix_mean_with_residual(double *out, double *residual,
                                                std::uint64_t count,
                                                std::size_t row_stride) const
{
    const std::size_t stride = row_stride ? row_stride : cols_;
    for (std::size_t j = 0; j < cols_; ++j) {
        for (std::size_t i = 0; i < rows_; ++i) {
            out[i * stride + j] =
                to_double_mean(&residual[i * stride + j], i, j, count);
        }
    }
}

bool
AccumulationMatrix::is_exactly_representable(std::size_t i, std::size_t j) const
{
    bool neg = false;
    const std::vector<limb_t> mag = magnitude(&neg, i, j);
    const std::size_t b = limbs::bit_length(mag.data(), mag.size());
    if (0 == b) return true;

    std::size_t ecol = 0, eslot = 0;
    locate(ecol, eslot, i, j);
    const long long exp = cols_state_[ecol].exponent;
    const std::size_t tz = trailing_zeros(mag);
    if (b - tz > 53) return false;                                 // too wide
    if (exp + static_cast<long long>(b) - 1 > 1023) return false;  // overflow
    if (exp + static_cast<long long>(tz) < -1074) return false;    // underflow
    return true;
}

std::string
AccumulationMatrix::to_exact_decimal(std::size_t i, std::size_t j) const
{
    bool neg = false;
    std::vector<limb_t> mag = magnitude(&neg, i, j);
    if (0 == limbs::bit_length(mag.data(), mag.size())) return "0";

    std::size_t ecol = 0, eslot = 0;
    locate(ecol, eslot, i, j);
    const long long exp = cols_state_[ecol].exponent;
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
AccumulationMatrix::reserve_column(std::size_t j, int exponent,
                                   std::size_t bits)
{
    if (j >= cols_)
        throw std::out_of_range("truesum: column index out of range");
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
AccumulationMatrix::reserve_for(const double *b, std::size_t count,
                                std::size_t row_stride)
{
    const std::size_t stride = row_stride ? row_stride : cols_;
    const std::size_t headroom = ceil_log2(std::max<std::size_t>(count, 1)) + 1;

    for (std::size_t j = 0; j < cols_; ++j) {
        bool any = false;
        long long low = 0, high = 0;
        const Column &cj = cols_state_[j];
        for (std::size_t k = 0; k < cj.rows; ++k) {
            const double v = b[(cj.first_row + k) * stride + j];
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
AccumulationMatrix::memory_bytes() const
{
    std::size_t total = sizeof(*this);
    for (const auto &c : cols_state_) {
        const std::size_t per_limb =
            padded_rows(c.rows) * sizeof(limb_t) + 512;
        total += sizeof(Column) + c.limbs.size() * per_limb +
                 c.bases.capacity() * sizeof(limb_t *);
    }
    return total;
}

std::string
AccumulationMatrix::describe() const
{
    std::ostringstream os;
    os << rows_ << " x " << cols_ << " column-block fixed-point matrix (";
    if (symmetric_) {
        os << "symmetric, " << (Uplo::Lower == uplo_ ? "lower" : "upper")
           << " stored, ";
    }
    os << memory_bytes() << " bytes, " << active_kernel() << " kernel)\n";
    for (std::size_t j = 0; j < cols_; ++j) {
        const Column &c = cols_state_[j];
        os << "  col " << j << ": exponent 2^" << c.exponent << ", width "
           << c.limbs.size() * kLimbBits << " bits";
        if (!c.initialized) os << " (unset)";
        os << "\n";
    }
    return os.str();
}

}  // namespace truesum
