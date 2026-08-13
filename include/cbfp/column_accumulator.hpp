// Exact accumulation of double-valued matrices into a column-blocked
// arbitrary-precision fixed-point matrix.
//
// Every column carries its own scale. Entry (i, j) is stored as a two's
// complement integer V of `column_limbs(j)` 64-bit limbs and denotes the exact
// rational value
//
//     A(i, j) = V * 2^column_exponent(j)
//
// Because a double is exactly m * 2^e for an integer m, every accumulation is
// exact provided the column exponent is low enough and the block is wide
// enough. Both adapt automatically: the exponent only ever decreases and the
// width only ever increases, so no intermediate precision is ever discarded.
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "cbfp/limb_column.hpp"
#include "cbfp/limbs.hpp"
#include "cbfp/survey.hpp"

namespace cbfp {

// Name of the kernel variant selected for this CPU ("scalar", "avx512").
const char *
active_kernel();

namespace detail {
class ThreadPool;
}

class ColumnBlockMatrix {
public:
    ColumnBlockMatrix(std::size_t rows, std::size_t cols);

    // Pre-sized from the surveys of every matrix that will be accumulated:
    // `surveys` is indexed by matrix and then by column, so its outer size is
    // how many will arrive and supplies the headroom term that would otherwise
    // come from counting them.
    //
    // Sized this way the accumulator never rescales, never widens, and never
    // surveys an incoming batch, because it already knows what is in one. The
    // metadata is taken on trust: a matrix reaching outside what was declared,
    // or more matrices than were described, voids the sizing and the sums are
    // then simply wrong. Only the count is checked, because that costs nothing.
    ColumnBlockMatrix(std::size_t rows, std::size_t cols,
                      const std::vector<std::vector<Survey>> &surveys);

    // Declared, not implicit: the worker pool is held by unique_ptr to an
    // incomplete type, so the destructor has to be defined where that type is.
    // The moves are spelled out because declaring a destructor would otherwise
    // suppress them.
    ~ColumnBlockMatrix();
    ColumnBlockMatrix(ColumnBlockMatrix &&) noexcept;
    ColumnBlockMatrix &operator=(ColumnBlockMatrix &&) noexcept;

    std::size_t rows() const { return rows_; }

    std::size_t cols() const { return cols_; }

    // --- exact accumulation ------------------------------------------------

    // A(i, j) += v. Exact. Throws std::domain_error on inf/NaN.
    void add(std::size_t i, std::size_t j, double v);

    // A(i, j) -= v. Exact.
    void sub(std::size_t i, std::size_t j, double v);

    // A += B, where B is a dense row-major rows() x cols() matrix whose rows
    // are `row_stride` doubles apart (0 means "tightly packed", i.e. cols()).
    void add_matrix(const double *b, std::size_t row_stride = 0);

    // A += scale * B, where scale is a power of two. Exact for any such scale
    // (including negative ones); this is the only scaling that stays exact
    // without widening the mantissa.
    void add_matrix_scaled_pow2(const double *b, int log2_scale,
                                std::size_t row_stride = 0);

    // A += B, where B is column-major: column j begins at b + j*col_stride and
    // its rows are contiguous (0 means tightly packed, i.e. rows()). This is
    // the layout the accumulator itself wants, so it avoids the staging copy
    // that a row-major input needs.
    void add_matrix_col_major(const double *b, std::size_t col_stride = 0);

    void add_matrix_col_major_scaled_pow2(const double *b, int log2_scale,
                                          std::size_t col_stride = 0);

    // A(., j) += v, where v is rows() contiguous doubles. Lets a caller stream
    // one column at a time, so only the accumulator has to be resident -- the
    // difference between holding a whole input matrix and holding 8*rows()
    // bytes of it.
    void add_column(std::size_t j, const double *v);

    void add_column_scaled_pow2(std::size_t j, const double *v, int log2_scale);

    void set_zero();

    // --- threading ---------------------------------------------------------

    // Spread accumulation across `n` worker threads, partitioned by column.
    // Columns are independent in storage, so this needs no locking and gives
    // bit-identical results -- exact accumulation is order independent, and
    // each column is touched by exactly one thread besides.
    //
    // 1 is the default and runs everything on the calling thread. Threads are
    // created here and kept, because creating them per call costs more than it
    // saves on small matrices.
    void set_threads(unsigned n);

    unsigned threads() const;

    // --- readback ----------------------------------------------------------

    // Correctly rounded (round-to-nearest, ties-to-even) double nearest to the
    // exact stored value. Overflows to +/-inf; underflows through the
    // subnormal range without double rounding.
    double to_double(std::size_t i, std::size_t j) const;

    // Fills a dense row-major rows() x cols() buffer with to_double() results.
    void to_matrix(double *out, std::size_t row_stride = 0) const;

    // The exact value as a decimal string, e.g. "0.1" accumulated once yields
    // 0.1000000000000000055511151231257827021181583404541015625. Never rounds.
    std::string to_exact_decimal(std::size_t i, std::size_t j) const;

    // True if the exact stored value is representable as a double with no
    // rounding at all.
    bool is_exactly_representable(std::size_t i, std::size_t j) const;

    bool is_zero(std::size_t i, std::size_t j) const;

    // Entry (i, j)'s stored two's complement limbs, little-endian, exactly as
    // held. Not a hot path: this exists so a second implementation can be
    // checked against this one. Comparing limbs is strictly stronger than
    // comparing to_double or to_exact_decimal, since identical limbs and an
    // identical column exponent imply identical everything downstream.
    std::vector<limbs::limb_t> entry_limbs(std::size_t i, std::size_t j) const;

    // --- column state / tuning ---------------------------------------------

    // Power-of-two weight of bit 0 of this column's stored integers.
    int column_exponent(std::size_t j) const { return cols_state_[j].exponent; }

    std::size_t column_bit_width(std::size_t j) const
    {
        return cols_state_[j].limbs.size() * limbs::kLimbBits;
    }

    std::size_t column_limbs(std::size_t j) const
    {
        return cols_state_[j].limbs.size();
    }

    // Pre-size a column so that accumulation performs no rescaling:
    // `exponent` is the lowest bit weight that will be needed, `bits` the
    // total width. Widening/lowering only; never discards precision.
    void reserve_column(std::size_t j, int exponent, std::size_t bits);

    // Pre-size every column from a matrix that is about to be accumulated
    // `count` times, so the accumulation itself never rescales. Purely an
    // optimization; results are identical without it.
    void reserve_for(const double *b, std::size_t count = 1,
                     std::size_t row_stride = 0);

    std::size_t memory_bytes() const;

    // Human-readable per-column exponent/width report.
    std::string describe() const;

private:
    struct Column {
        int exponent = 0;
        bool initialized = false;  // false until the first add fixes the scale

        // Width is derived rather than measured: no entry can exceed
        // `count * 2^max_addend_bits`, so that bound plus a sign bit says how
        // many limbs are needed, and the kernels can then skip overflow checks.
        std::size_t max_addend_bits = 0;
        std::size_t add_count = 0;

        std::vector<LimbColumn> limbs;       // limbs[k] = k-th limb of a row
        std::vector<limbs::limb_t *> bases;  // limbs[k].data(), for kernels
    };

    // Column j starts at b + j*column_step, with element i at + i*row_step.
    // Row-major is (1, row_stride); column-major is (col_stride, 1).
    void accumulate_columns(const double *b, std::size_t column_step,
                            std::size_t row_step, int log2_scale);

    // The columns in [begin, end) of one such matrix, staged through the
    // buffer belonging to worker `slot`.
    void accumulate_column_range(const double *b, std::size_t column_step,
                                 std::size_t row_step, int log2_scale,
                                 std::size_t begin, std::size_t end,
                                 unsigned slot);

    // Accumulates one column that is already contiguous.
    void accumulate_column(std::size_t j, const double *column, int log2_scale);

    void check_index(std::size_t i, std::size_t j) const;
    void report_contradiction(std::size_t j, unsigned flags) const;
    void rebuild_bases(Column &c);
    void ensure_limb_count(Column &c, std::size_t limbs_needed);
    void fit_column(Column &c);
    void rescale(Column &c, int new_exponent);
    bool column_is_zero(const Column &c) const;
    void accumulate(std::size_t i, std::size_t j, double v, bool negate,
                    int log2_scale);

    // Absolute value of entry (i, j) as a magnitude limb vector, plus its
    // sign. Gathers across limb positions, so this is a readback path only.
    std::vector<limbs::limb_t> magnitude(std::size_t i, std::size_t j,
                                         bool *negative) const;

    std::size_t rows_;
    std::size_t cols_;
    std::vector<Column> cols_state_;

    // Set when the constructor was given surveys. The accumulation path then
    // skips the survey, the rescale, the widen and the running width
    // bookkeeping -- all of which exist to discover what it was already told.
    bool presized_ = false;
    std::size_t declared_matrices_ = 0;
    std::size_t submitted_matrices_ = 0;

    // A column of a row-major matrix is strided, and a strided vector gather
    // costs more than the decomposition it feeds. Staging the column here once
    // lets both passes read contiguously. One buffer per worker, since that
    // was the only state the columns shared.
    std::vector<std::vector<double>> column_buffers_;

    // Workers, parked between calls. Held by pointer so this header does not
    // drag in <thread> and friends.
    std::unique_ptr<detail::ThreadPool> pool_;
};

// Decomposes a finite double into an exact odd mantissa and exponent:
// v == mantissa * 2^exponent, with mantissa odd (or zero, when v is +/-0).
// For +/-0 the mantissa is 0 and the exponent carries no meaning; `negative`
// still reflects the sign bit, so -0.0 reports negative.
struct DoubleParts {
    std::uint64_t mantissa;  // magnitude, at most 53 significant bits, odd
    bool negative;
    int exponent;
};

DoubleParts
decompose(double v);

}  // namespace cbfp
