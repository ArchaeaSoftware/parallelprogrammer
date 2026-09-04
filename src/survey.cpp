#include "truesum/survey.hpp"

#include <climits>
#include <vector>

#include "kernels.hpp"

namespace truesum {

Survey
survey_column(const double *values, std::size_t rows)
{
    // A producer has no column exponent to compare against, so the floor is
    // one nothing clears and the exact true-ulp minimum is always computed.
    return kernels::survey()(values, rows, LLONG_MAX);
}

void
survey_matrix_col_major(Survey *out, const double *b, std::size_t rows,
                        std::size_t cols, std::size_t col_stride)
{
    const std::size_t stride = col_stride ? col_stride : rows;
    for (std::size_t j = 0; j < cols; ++j) {
        out[j] = survey_column(b + j * stride, rows);
    }
}

void
survey_matrix(Survey *out, const double *b, std::size_t rows, std::size_t cols,
              std::size_t row_stride)
{
    const std::size_t stride = row_stride ? row_stride : cols;
    // The kernels want a contiguous column and a strided gather costs more
    // than the work it feeds, so each column is staged once -- the same
    // trade the accumulator makes for row-major input.
    std::vector<double> column(rows);
    for (std::size_t j = 0; j < cols; ++j) {
        for (std::size_t i = 0; i < rows; ++i) column[i] = b[i * stride + j];
        out[j] = survey_column(column.data(), rows);
    }
}

}  // namespace truesum
