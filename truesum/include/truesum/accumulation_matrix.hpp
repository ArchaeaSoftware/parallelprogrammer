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

#include "truesum/limb_column.hpp"
#include "truesum/limbs.hpp"
#include "truesum/survey.hpp"
#include "truesum/uplo.hpp"

namespace truesum {

// Name of the kernel variant selected for this CPU ("scalar", "avx512").
const char *
active_kernel();

namespace detail {
class ThreadPool;
}

class AccumulationMatrix {
public:
    AccumulationMatrix(std::size_t rows, std::size_t cols);

    // Symmetric n x n, storing one triangle. A(i, j) and A(j, i) are the same
    // stored entry rather than two that happen to agree, which matters here
    // beyond halving the memory: every column carries its own exponent and
    // width, so in full storage the two halves would hold equal values in
    // different limbs, and symmetry would be a property the data had to keep
    // earning. Stored once, it cannot drift.
    //
    // Accumulating a matrix reads only the stored triangle of the input, so
    // the input traffic halves as well. The input is taken to be symmetric;
    // that is not checked, because checking it means reading the half this
    // exists to avoid reading.
    AccumulationMatrix(std::size_t n, Uplo uplo);

    // Pre-sized from the surveys of every matrix that will be accumulated:
    // `surveys` is indexed by matrix and then by column, so its outer size is
    // how many will arrive and supplies the headroom term that would otherwise
    // come from counting them.
    //
    // Sized this way the accumulation matrix never rescales, never widens, and
    // never surveys an incoming batch, because it already knows what is in one.
    // The metadata is taken on trust: a matrix reaching outside what was
    // declared, or more matrices than were described, voids the sizing and the
    // sums are then simply wrong. Only the count is checked, because that costs
    // nothing.
    AccumulationMatrix(std::size_t rows, std::size_t cols,
                       const std::vector<std::vector<Survey>> &surveys);

    // Both at once. Each survey describes the stored triangle of one matrix,
    // column by column: entry j must cover the column_rows(j) values beginning
    // at row column_first_row(j), since those are the only ones that will be
    // read.
    AccumulationMatrix(std::size_t n, Uplo uplo,
                       const std::vector<std::vector<Survey>> &surveys);

    // Declared, not implicit: the worker pool is held by unique_ptr to an
    // incomplete type, so the destructor has to be defined where that type is.
    // The moves are spelled out because declaring a destructor would otherwise
    // suppress them.
    ~AccumulationMatrix();
    AccumulationMatrix(AccumulationMatrix &&) noexcept;
    AccumulationMatrix &operator=(AccumulationMatrix &&) noexcept;

    std::size_t rows() const { return rows_; }

    std::size_t cols() const { return cols_; }

    bool symmetric() const { return symmetric_; }

    // Only meaningful when symmetric().
    Uplo uplo() const { return uplo_; }

    // How many entries column j actually stores: rows() unless symmetric, and
    // then n-j for Lower or j+1 for Upper. This is the length add_column()
    // wants, and the length of the slice of an input column that is read.
    std::size_t column_rows(std::size_t j) const;

    // The logical row index that column j's first stored entry corresponds to:
    // 0 unless symmetric and Lower, where it is j.
    std::size_t column_first_row(std::size_t j) const;

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
    // the layout the accumulation matrix itself wants, so it avoids the staging
    // copy that a row-major input needs.
    void add_matrix_col_major(const double *b, std::size_t col_stride = 0);

    void add_matrix_col_major_scaled_pow2(const double *b, int log2_scale,
                                          std::size_t col_stride = 0);

    // A += B[0] + B[1] + ... + B[count-1], all in one pass over the
    // accumulation matrix. `b` is an array of `count` pointers, each to a
    // column-major matrix laid out as add_matrix_col_major expects.
    //
    // Identical results to calling add_matrix_col_major once per matrix --
    // exact accumulation does not care about order or grouping. What changes
    // is traffic. One at a time, each batch reads and writes every limb it
    // touches; folded, the limbs are read once, all `count` addends applied in
    // registers, and written once. That is a saving in accumulation matrix
    // traffic only, so it shows up where accumulation matrix traffic is the
    // bound.
    //
    // Which `count` to pass is not obvious, because two costs pull the other
    // way: `count` input columns are `count` concurrent streams rather than
    // one, and a row's limbs have to stay live across the whole fold. Measured
    // at 64 columns, Gelem/s, medians of three:
    //
    //                        seq   K=2   K=4   K=8
    //   65536 rows, 8 thr   1.14  1.84  2.65  3.01   <- DRAM-bound
    //   any shape, 1 thr    0.57  0.69  0.59  0.49
    //
    // So: past L3 and threaded, fold as much as you have -- 2.6x at K=8, and
    // still rising. Single-threaded, K=2 is worth 1.2x and beyond that the
    // extra streams cost more than the traffic saved. Inside L3 on many
    // threads the run-to-run spread swamps the difference, and the ordinary
    // entry point is the simpler one.
    //
    // Column-major only, and deliberately: folding a row-major batch would
    // mean staging `count` columns at once per worker, and writing and
    // re-reading that staging is the traffic this exists to avoid.
    //
    // Columns wider than the fold's register budget fall back to one batch at
    // a time, which costs nothing: the single-batch kernel stops as soon as a
    // carry dies, while the fold must write back every limb it loaded.
    void add_matrices_col_major(const double *const *b, std::size_t count,
                                std::size_t col_stride = 0);

    // A(., j) += v, where v is column_rows(j) contiguous doubles. Lets a caller
    // stream one column at a time, so only the accumulation matrix has to be
    // resident
    // -- the difference between holding a whole input matrix and holding
    // 8*rows() bytes of it.
    //
    // When symmetric, v is the stored slice of the column, not the whole
    // column: element k of v is logical entry (column_first_row(j) + k, j).
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
    // denormal range without double rounding.
    double to_double(std::size_t i, std::size_t j) const;

    // The same, and the residual through the passback, which leads because the
    // function writes through it. Two entry points rather than one defaulted
    // parameter: a caller who wants the residual says so by which one they
    // call, and a caller who does not is not made to write nullptr.
    //
    // Writes what the rounding discarded:
    // exactly (stored value - returned double), itself correctly rounded to a
    // double. The subtraction happens in the accumulation matrix's own fixed
    // point, so the residual is the true difference and not an estimate --
    // forming it in floating point would be the cancellation this container
    // exists to avoid.
    //
    // Properties worth relying on:
    //   - the returned double is unaffected by asking for the residual;
    //   - |residual| <= half an ulp of the returned double, and the residual
    //     is negative exactly when the rounding went up, so the two together
    //     are a non-overlapping two-term expansion carrying about 106 bits;
    //   - the residual is +0.0 when the value is exactly representable, which
    //     is_exactly_representable() reports independently;
    //   - a value outside double's range returns +/-inf with a residual of
    //     0.0, the difference being unrepresentable.
    double to_double(double *residual, std::size_t i, std::size_t j) const;

    // Fills a dense row-major rows() x cols() buffer with to_double() results.
    void to_matrix(double *out, std::size_t row_stride = 0) const;

    // to_matrix, plus the residual of every entry written to a second buffer
    // of the same shape and stride. residual[i*stride + j] belongs to
    // out[i*stride + j].
    void to_matrix_with_residual(double *out, double *residual,
                                 std::size_t row_stride = 0) const;

    // --- means -------------------------------------------------------------

    // The stored sum divided by `count`, correctly rounded: the double nearest
    // the exact rational (stored value / count), round-to-nearest ties-to-
    // even, over the full range and through the denormals, exactly as
    // to_double is for the sum itself. A power-of-two count is exact scaling;
    // any other count is where a naive `to_double(i, j) / count` rounds twice,
    // and this does not.
    //
    // The count is the caller's, not the accumulation matrix's. The
    // accumulation matrix does not know how many values a cell's sum stands
    // for: add_matrix, add_column and add all feed the same cell, a scaled add
    // is a weighted one, and the internal add counter is a width bound that
    // resets on rescale. Passing the count also makes this a general exact
    // division by an integer, which is what a weighted mean needs.
    //
    // A count of zero throws std::domain_error.
    double to_double_mean(std::size_t i, std::size_t j,
                          std::uint64_t count) const;

    // The same, with the residual: exactly (stored value / count - returned
    // double), correctly rounded. The quotient is not a fixed-point value, so
    // the difference is formed as an integer numerator over the same count
    // and divided once more; nothing is rounded before the final step. The
    // properties listed for to_double's residual all hold here, with
    // is_exactly_representable read as "the residual is +0.0".
    double to_double_mean(double *residual, std::size_t i, std::size_t j,
                          std::uint64_t count) const;

    // to_matrix and to_matrix_with_residual, over the mean.
    void to_matrix_mean(double *out, std::uint64_t count,
                        std::size_t row_stride = 0) const;

    void to_matrix_mean_with_residual(double *out, double *residual,
                                      std::uint64_t count,
                                      std::size_t row_stride = 0) const;

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
        // many limbs are needed, and the carry chain then needs no overflow
        // test below its top limb.
        std::size_t max_addend_bits = 0;
        std::size_t add_count = 0;

        // How many entries this column stores, and which logical row the first
        // of them is. Held per column rather than derived from uplo because a
        // triangular column has its own length, every allocation, rescale and
        // sweep needs it, and the accumulate loop should not branch on shape.
        std::size_t rows = 0;
        std::size_t first_row = 0;

        std::vector<LimbColumn> limbs;       // limbs[k] = k-th limb of a row
        std::vector<limbs::limb_t *> bases;  // limbs[k].data(), for kernels
    };

    // Logical (i, j) to the entry that actually holds it: `col` is the storing
    // column and `slot` the index within it. The identity unless symmetric,
    // where the triangle folds one index onto the other.
    void locate(std::size_t &col, std::size_t &slot, std::size_t i,
                std::size_t j) const;

    // Shared tail of the two survey-taking constructors.
    void reserve_from_surveys(const std::vector<std::vector<Survey>> &surveys);

    // The columns belonging to worker `slot`. Even by column count normally;
    // by stored entries when symmetric, since a triangular column's work runs
    // from n down to 1 and splitting on count alone would leave one worker
    // with most of the matrix.
    void partition_columns(std::size_t &begin, std::size_t &end,
                           unsigned slot) const;

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

    // The fold's counterparts: one column of every matrix at once.
    void fold_column_range(const double *const *b, std::size_t count,
                           std::size_t col_stride, std::size_t begin,
                           std::size_t end);
    void fold_column(std::size_t j, const double *const *columns,
                     std::size_t count);

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
    std::vector<limbs::limb_t> magnitude(bool *negative, std::size_t i,
                                         std::size_t j) const;

    std::size_t rows_;
    std::size_t cols_;
    bool symmetric_ = false;
    Uplo uplo_ = Uplo::Lower;
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

}  // namespace truesum
