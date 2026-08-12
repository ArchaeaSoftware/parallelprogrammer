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
#include <string>
#include <vector>

#include "cbfp/limbs.hpp"

namespace cbfp {

class ColumnBlockMatrix {
public:
    ColumnBlockMatrix(std::size_t rows, std::size_t cols);

    std::size_t rows() const { return rows_; }

    std::size_t cols() const { return cols_; }

    // --- exact accumulation ------------------------------------------------

    // A(i, j) += v. Exact. Throws std::domain_error on inf/NaN.
    void add(std::size_t i, std::size_t j, double v);

    // A(i, j) -= v. Exact.
    void sub(std::size_t i, std::size_t j, double v);

    // A += B, where B is a dense row-major rows() x cols() matrix whose rows
    // are `row_stride` doubles apart (0 means "tightly packed", i.e. cols()).
    void add_matrix(const double* b, std::size_t row_stride = 0);

    // A += scale * B, where scale is a power of two. Exact for any such scale
    // (including negative ones); this is the only scaling that stays exact
    // without widening the mantissa.
    void add_matrix_scaled_pow2(const double* b, int log2_scale,
                                std::size_t row_stride = 0);

    void set_zero();

    // --- readback ----------------------------------------------------------

    // Correctly rounded (round-to-nearest, ties-to-even) double nearest to the
    // exact stored value. Overflows to +/-inf; underflows through the
    // subnormal range without double rounding.
    double to_double(std::size_t i, std::size_t j) const;

    // Fills a dense row-major rows() x cols() buffer with to_double() results.
    void to_matrix(double* out, std::size_t row_stride = 0) const;

    // The exact value as a decimal string, e.g. "0.1" accumulated once yields
    // 0.1000000000000000055511151231257827021181583404541015625. Never rounds.
    std::string to_exact_decimal(std::size_t i, std::size_t j) const;

    // True if the exact stored value is representable as a double with no
    // rounding at all.
    bool is_exactly_representable(std::size_t i, std::size_t j) const;

    bool is_zero(std::size_t i, std::size_t j) const;

    // --- column state / tuning ---------------------------------------------

    // Power-of-two weight of bit 0 of this column's stored integers.
    int column_exponent(std::size_t j) const { return cols_state_[j].exponent; }

    std::size_t column_bit_width(std::size_t j) const
    {
        return cols_state_[j].limbs * limbs::kLimbBits;
    }

    std::size_t column_limbs(std::size_t j) const
    {
        return cols_state_[j].limbs;
    }

    // Pre-size a column so that accumulation performs no rescaling:
    // `exponent` is the lowest bit weight that will be needed, `bits` the
    // total width. Widening/lowering only; never discards precision.
    void reserve_column(std::size_t j, int exponent, std::size_t bits);

    // Pre-size every column from a matrix that is about to be accumulated
    // `count` times, so the accumulation itself never rescales. Purely an
    // optimization; results are identical without it.
    void reserve_for(const double* b, std::size_t count = 1,
                     std::size_t row_stride = 0);

    std::size_t memory_bytes() const;

    // Human-readable per-column exponent/width report.
    std::string describe() const;

private:
    struct Column {
        int exponent = 0;
        std::size_t limbs = 1;
        bool initialized = false;  // false until the first add fixes the scale
        std::vector<limbs::limb_t> data;
    };

    void check_index(std::size_t i, std::size_t j) const;
    void ensure_limbs(Column& c, std::size_t limbs_needed);
    void rescale(Column& c, int new_exponent);
    void accumulate(std::size_t i, std::size_t j, double v, bool negate,
                    int log2_scale);

    limbs::limb_t* entry(Column& c, std::size_t i)
    {
        return c.data.data() + i * c.limbs;
    }

    const limbs::limb_t* entry(const Column& c, std::size_t i) const
    {
        return c.data.data() + i * c.limbs;
    }

    // Absolute value of entry (i, j) as a magnitude limb vector, plus its
    // sign.
    std::vector<limbs::limb_t> magnitude(std::size_t i, std::size_t j,
                                         bool* negative) const;

    std::size_t rows_;
    std::size_t cols_;
    std::vector<Column> cols_state_;
};

// Decomposes a finite double into an exact odd mantissa and exponent:
// v == mantissa * 2^exponent, with mantissa odd (or zero, when v is +/-0).
struct DoubleParts {
    std::uint64_t mantissa;  // magnitude, at most 53 significant bits, odd
    bool negative;
    int exponent;
};

DoubleParts decompose(double v);

}  // namespace cbfp
