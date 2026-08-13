// What an accumulator needs to know about a matrix before it can size itself
// for one.
//
// This is deliberately public, and small. Sizing a column requires the lowest
// bit weight its values will reach and the highest they will occupy, which
// means reading every element -- work that whoever produced the matrix has
// already done, with the values in registers. A survey is four numbers per
// column against eight bytes per element: 1 KB beside a 33.6 MB matrix at
// 65536x64, one part in 32768. So it can travel with the matrix rather than
// being recomputed from it, and an accumulator given the surveys of every
// matrix it will receive can be allocated once, exactly, and never rescale or
// widen afterwards.
#pragma once

#include <cstddef>

namespace cbfp {

struct Survey {
    // Lowest true-ulp exponent among the nonzero values -- the exponent after
    // normalizing each mantissa to odd, not the frexp exponent, so a column of
    // integers reports 0 rather than -52.
    long long min_exponent;

    // One past the highest bit position any value occupies, as a power of two.
    long long max_top;

    bool any;        // false if every value was zero
    bool nonfinite;  // true if any value was inf or NaN
};

// Surveys one column: `rows` contiguous doubles.
Survey
survey_column(const double *values, std::size_t rows);

// Surveys every column of a column-major matrix, writing `cols` entries.
// Column j begins at b + j*col_stride and its rows are contiguous; 0 means
// tightly packed.
void
survey_matrix_col_major(const double *b, std::size_t rows, std::size_t cols,
                        Survey *out, std::size_t col_stride = 0);

// The same for a row-major matrix, whose rows are `row_stride` doubles apart
// (0 means tightly packed, i.e. cols). Each column is strided, so this is the
// slower form; it exists so a producer holding a row-major matrix need not
// transpose one just to describe it.
void
survey_matrix(const double *b, std::size_t rows, std::size_t cols, Survey *out,
              std::size_t row_stride = 0);

}  // namespace cbfp
