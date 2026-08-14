#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "cbfp/column_accumulator.hpp"
#include "cbfp/survey.hpp"

namespace {

int g_failures = 0;
int g_checks = 0;

void
check(bool ok, const char *expr, const char *file, int line)
{
    ++g_checks;
    if (!ok) {
        ++g_failures;
        std::printf("FAIL %s:%d: %s\n", file, line, expr);
    }
}

#define CHECK(expr) check((expr), #expr, __FILE__, __LINE__)

void
check_eq_double(double got, double want, const char *what, int line)
{
    ++g_checks;
    const bool ok = (std::isnan(got) && std::isnan(want)) ||
                    (got == want && std::signbit(got) == std::signbit(want));
    if (!ok) {
        ++g_failures;
        std::printf("FAIL %s:%d: %s: got %.17g, want %.17g\n", __FILE__, line,
                    what, got, want);
    }
}

#define CHECK_DOUBLE(got, want) check_eq_double((got), (want), #got, __LINE__)

void
check_eq_str(const std::string &got, const std::string &want, const char *what,
             int line)
{
    ++g_checks;
    if (got != want) {
        ++g_failures;
        std::printf("FAIL %s:%d: %s:\n  got  %s\n  want %s\n", __FILE__, line,
                    what, got.c_str(), want.c_str());
    }
}

#define CHECK_STR(got, want) check_eq_str((got), (want), #got, __LINE__)

// ---------------------------------------------------------------------------

void
test_decompose()
{
    auto p = cbfp::decompose(1.0);
    CHECK(p.mantissa == 1 && p.exponent == 0 && !p.negative);

    p = cbfp::decompose(-8.0);
    CHECK(p.mantissa == 1 && p.exponent == 3 && p.negative);

    // Both zeros give a zero mantissa; the sign bit still comes through, but
    // carries no meaning because accumulate() returns early on a zero mantissa.
    p = cbfp::decompose(0.0);
    CHECK(p.mantissa == 0 && !p.negative);

    p = cbfp::decompose(-0.0);
    CHECK(p.mantissa == 0 && p.negative);

    p = cbfp::decompose(0.5);
    CHECK(p.mantissa == 1 && p.exponent == -1);

    // Smallest subnormal is exactly 2^-1074.
    p = cbfp::decompose(std::numeric_limits<double>::denorm_min());
    CHECK(p.mantissa == 1 && p.exponent == -1074);

    // Every mantissa comes back odd, and the value reconstructs exactly.
    std::mt19937_64 rng(1234);
    for (int n = 0; n < 2000; ++n) {
        double v;
        std::uint64_t bits = rng();
        std::memcpy(&v, &bits, sizeof v);
        if (!std::isfinite(v)) continue;

        p = cbfp::decompose(v);
        if (p.mantissa == 0) {
            CHECK(v == 0.0);
            continue;
        }
        CHECK((p.mantissa & 1) == 1);
        const double rebuilt =
            std::ldexp(static_cast<double>(p.mantissa), p.exponent);
        CHECK_DOUBLE(p.negative ? -rebuilt : rebuilt, v);
    }
}

void
test_single_value_roundtrip()
{
    const double values[] = {
        0.0,
        -0.0,
        1.0,
        -1.0,
        0.1,
        -0.1,
        3.14159265358979,
        std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max(),
        std::numeric_limits<double>::min(),         // smallest normal
        std::numeric_limits<double>::denorm_min(),  // smallest subnormal
        -std::numeric_limits<double>::denorm_min(),
        std::ldexp(1.0, 1000),
        std::ldexp(1.0, -1000),
        std::ldexp(4503599627370495.0, -1074),  // largest subnormal
    };

    for (double v : values) {
        cbfp::ColumnBlockMatrix m(1, 1);
        m.add(0, 0, v);
        // -0.0 accumulates as an exact zero; the sign of zero is not tracked.
        const double want = (v == 0.0) ? 0.0 : v;
        CHECK_DOUBLE(m.to_double(0, 0), want);
        CHECK(m.is_exactly_representable(0, 0));
    }
}

void
test_catastrophic_cancellation()
{
    // The classic case: naive double summation returns 0, the exact answer
    // is 1.
    const double terms[] = {1e300, 1.0, -1e300};
    double naive = 0.0;
    for (double t : terms) naive += t;
    CHECK_DOUBLE(naive, 0.0);

    cbfp::ColumnBlockMatrix m(1, 1);
    for (double t : terms) m.add(0, 0, t);
    CHECK_DOUBLE(m.to_double(0, 0), 1.0);
    CHECK_STR(m.to_exact_decimal(0, 0), "1");
}

void
test_repeated_tenth()
{
    // double(0.1) is exactly 3602879701896397 * 2^-55, so ten of them sum to
    // exactly 1 + 2^-54 -- a quarter ulp above 1, which rounds back to 1.0.
    // Naive summation instead drifts one ulp below.
    cbfp::ColumnBlockMatrix m(1, 1);
    double naive = 0.0;
    for (int i = 0; i < 10; ++i) {
        m.add(0, 0, 0.1);
        naive += 0.1;
    }
    CHECK_DOUBLE(naive, 0.99999999999999989);
    CHECK_DOUBLE(m.to_double(0, 0), 1.0);
    CHECK_STR(m.to_exact_decimal(0, 0),
              "1.000000000000000055511151231257827021181583404541015625");
}

void
test_exact_decimal()
{
    cbfp::ColumnBlockMatrix m(1, 4);
    m.add(0, 0, 0.1);
    CHECK_STR(m.to_exact_decimal(0, 0),
              "0.1000000000000000055511151231257827021181583404541015625");

    m.add(0, 1, -0.3);
    CHECK_STR(m.to_exact_decimal(0, 1),
              "-0.299999999999999988897769753748434595763683319091796875");

    m.add(0, 2, 1024.0);
    CHECK_STR(m.to_exact_decimal(0, 2), "1024");

    m.add(0, 3, std::numeric_limits<double>::denorm_min());
    const std::string tiny = m.to_exact_decimal(0, 3);
    CHECK(tiny.size() == 1076);  // "0." + 1074 digits
    CHECK(tiny.compare(0, 2, "0.") == 0);
    CHECK(tiny.back() == '5');  // 2^-1074 ends in ...625
}

void
test_exponent_and_width_tracking()
{
    cbfp::ColumnBlockMatrix m(1, 1);
    m.add(0, 0, 1.0);
    CHECK(m.column_exponent(0) == 0);

    // A much smaller value forces the column scale down; nothing is lost.
    m.add(0, 0, std::ldexp(1.0, -200));
    CHECK(m.column_exponent(0) == -200);
    CHECK(m.column_bit_width(0) >= 201);

    // A much larger value forces the block wider.
    m.add(0, 0, std::ldexp(1.0, 200));
    CHECK(m.column_exponent(0) == -200);
    CHECK(m.column_bit_width(0) >= 401);

    m.sub(0, 0, std::ldexp(1.0, 200));
    m.sub(0, 0, std::ldexp(1.0, -200));
    CHECK_DOUBLE(m.to_double(0, 0), 1.0);

    m.sub(0, 0, 1.0);
    CHECK(m.is_zero(0, 0));
    CHECK_DOUBLE(m.to_double(0, 0), 0.0);
}

void
test_full_double_range()
{
    // Span the entire binade range in a single column: 2^-1074 up to 2^1023.
    cbfp::ColumnBlockMatrix m(1, 1);
    m.add(0, 0, std::ldexp(1.0, 1023));
    m.add(0, 0, std::numeric_limits<double>::denorm_min());
    CHECK(m.column_bit_width(0) >= 1024 + 1074);

    // The tiny term is far below the ulp of the large one, so the rounded
    // readback is the large term, but the stored value still knows about it.
    CHECK_DOUBLE(m.to_double(0, 0), std::ldexp(1.0, 1023));
    CHECK(!m.is_exactly_representable(0, 0));

    m.sub(0, 0, std::ldexp(1.0, 1023));
    CHECK_DOUBLE(m.to_double(0, 0), std::numeric_limits<double>::denorm_min());
    CHECK(m.is_exactly_representable(0, 0));
}

void
test_rounding_ties_to_even()
{
    {  // 1 + 2^-53 is exactly halfway between 1 and nextafter(1); ties to
       // even.
        cbfp::ColumnBlockMatrix m(1, 1);
        m.add(0, 0, 1.0);
        m.add(0, 0, std::ldexp(1.0, -53));
        CHECK_DOUBLE(m.to_double(0, 0), 1.0);
    }
    {  // One bit above the tie rounds up.
        cbfp::ColumnBlockMatrix m(1, 1);
        m.add(0, 0, 1.0);
        m.add(0, 0, std::ldexp(1.0, -53));
        m.add(0, 0, std::ldexp(1.0, -105));
        CHECK_DOUBLE(m.to_double(0, 0), std::nextafter(1.0, 2.0));
    }
    {  // Tie with an odd mantissa rounds up (to even).
        cbfp::ColumnBlockMatrix m(1, 1);
        const double odd = std::nextafter(1.0, 2.0);  // 1 + 2^-52
        m.add(0, 0, odd);
        m.add(0, 0, std::ldexp(1.0, -53));
        CHECK_DOUBLE(m.to_double(0, 0), std::nextafter(odd, 2.0));
    }
    {  // 2^-1075 is halfway between 0 and the smallest subnormal: ties to
       // zero. It is below anything a double can hold, so it is reached by
       // accumulating the smallest subnormal with an exact power-of-two scale.
        cbfp::ColumnBlockMatrix m(1, 1);
        double sub = std::numeric_limits<double>::denorm_min();  // 2^-1074
        m.add_matrix_scaled_pow2(&sub, -1);                      // += 2^-1075
        CHECK_DOUBLE(m.to_double(0, 0), 0.0);
        CHECK(!m.is_zero(0, 0));  // exact value kept, only readback rounds
        CHECK(m.column_exponent(0) == -1075);

        // Just past that tie rounds up to the smallest subnormal.
        m.add_matrix_scaled_pow2(&sub, -2);  // += 2^-1076
        CHECK_DOUBLE(m.to_double(0, 0),
                     std::numeric_limits<double>::denorm_min());
    }
    {  // Two subnormals that sum into the smallest normal.
        cbfp::ColumnBlockMatrix m(1, 1);
        const double largest_sub = std::ldexp(4503599627370495.0, -1074);
        m.add(0, 0, largest_sub);
        m.add(0, 0, std::numeric_limits<double>::denorm_min());
        CHECK_DOUBLE(m.to_double(0, 0), std::numeric_limits<double>::min());
        CHECK(m.is_exactly_representable(0, 0));
    }
    {  // Overflow of the double range on readback only.
        cbfp::ColumnBlockMatrix m(1, 1);
        const double big = std::numeric_limits<double>::max();
        m.add(0, 0, big);
        m.add(0, 0, big);
        CHECK_DOUBLE(m.to_double(0, 0),
                     std::numeric_limits<double>::infinity());
        m.sub(0, 0, big);
        CHECK_DOUBLE(m.to_double(0, 0), big);  // and it comes back exactly
    }
}

void
test_add_then_subtract_is_zero()
{
    // Order-independence and exactness: adding a random set and then
    // subtracting the same values in a different order must land on exact
    // zero.
    std::mt19937_64 rng(20250812);
    std::uniform_int_distribution<int> exp_dist(-300, 300);
    std::uniform_real_distribution<double> man_dist(-1.0, 1.0);

    const std::size_t rows = 7, cols = 5;
    cbfp::ColumnBlockMatrix m(rows, cols);

    std::vector<std::vector<double>> per_cell(rows * cols);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            for (int n = 0; n < 40; ++n) {
                const double v = std::ldexp(man_dist(rng), exp_dist(rng));
                per_cell[i * cols + j].push_back(v);
                m.add(i, j, v);
            }
        }
    }
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            auto &cell = per_cell[i * cols + j];
            std::shuffle(cell.begin(), cell.end(), rng);
            for (double v : cell) m.sub(i, j, v);
            CHECK(m.is_zero(i, j));
        }
    }
}

void
test_matrix_accumulation()
{
    const std::size_t rows = 3, cols = 4;
    // clang-format off
    std::vector<double> b = {
        1e300,  1.0,  0.1, -7.0,
        1.0,    1e-8, 0.1,  2.5,
        -1e300, 3.0,  0.1,  0.0,
    };
    // clang-format on

    cbfp::ColumnBlockMatrix m(rows, cols);
    m.add_matrix(b.data());
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            CHECK_DOUBLE(m.to_double(i, j), b[i * cols + j]);
        }
    }

    // Accumulate the same matrix 1000 more times; every cell is exactly
    // 1001x.
    for (int n = 0; n < 1000; ++n) m.add_matrix(b.data());
    std::vector<double> out(rows * cols);
    m.to_matrix(out.data());
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            const double want = 1001.0 * b[i * cols + j];
            CHECK_DOUBLE(out[i * cols + j], want);
        }
    }

    // Power-of-two scaling is exact.
    cbfp::ColumnBlockMatrix h(rows, cols);
    h.add_matrix_scaled_pow2(b.data(), -3);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            CHECK_DOUBLE(h.to_double(i, j), b[i * cols + j] / 8.0);
        }
    }
}

void
test_reserve_avoids_rescaling()
{
    const std::size_t rows = 4, cols = 3;
    // clang-format off
    std::vector<double> b = {
        1.0,  1e-30, 1e20,
        2.0,  2e-30, 2e20,
        0.5,  4e-30, 4e20,
        0.25, 8e-30, 8e20,
    };
    // clang-format on

    cbfp::ColumnBlockMatrix lazy(rows, cols);
    cbfp::ColumnBlockMatrix eager(rows, cols);
    eager.reserve_for(b.data(), 500);

    for (int n = 0; n < 500; ++n) {
        lazy.add_matrix(b.data());
        eager.add_matrix(b.data());
    }
    for (std::size_t j = 0; j < cols; ++j) {
        CHECK(lazy.column_exponent(j) == eager.column_exponent(j));
    }
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            CHECK_STR(lazy.to_exact_decimal(i, j),
                      eager.to_exact_decimal(i, j));
        }
    }
}

void
test_column_independence()
{
    // Column 0 gets a huge dynamic range, column 1 stays cheap. The whole
    // point of per-column scaling is that column 1 does not pay for column 0.
    cbfp::ColumnBlockMatrix m(2, 2);
    m.add(0, 0, std::ldexp(1.0, 900));
    m.add(0, 0, std::ldexp(1.0, -900));
    m.add(0, 1, 1.0);
    m.add(1, 1, 2.0);

    CHECK(m.column_bit_width(0) >= 1800);
    CHECK(m.column_bit_width(1) <= 128);
    CHECK(m.column_exponent(1) == 0);
}

void
test_errors()
{
    cbfp::ColumnBlockMatrix m(2, 2);
    bool threw = false;
    try {
        m.add(0, 0, std::numeric_limits<double>::infinity());
    } catch (const std::domain_error &) {
        threw = true;
    }
    CHECK(threw);

    threw = false;
    try {
        m.add(5, 0, 1.0);
    } catch (const std::out_of_range &) {
        threw = true;
    }
    CHECK(threw);
}

void
test_many_small_into_large()
{
    // A million values whose ulp is far below the running total's ulp. Naive
    // summation stalls completely; this does not.
    cbfp::ColumnBlockMatrix m(1, 1);
    const double big = std::ldexp(1.0, 60);
    const double small = 1.0;
    m.add(0, 0, big);
    double naive = big;
    for (int n = 0; n < 1000000; ++n) {
        m.add(0, 0, small);
        naive += small;
    }
    CHECK_DOUBLE(naive, big);  // every increment is swallowed by rounding
    CHECK_DOUBLE(m.to_double(0, 0), big + 1000000.0);
    m.sub(0, 0, big);
    CHECK_DOUBLE(m.to_double(0, 0), 1000000.0);
}

void
test_reuse_after_zero()
{
    cbfp::ColumnBlockMatrix m(1, 1);
    m.add(0, 0, 1.0);
    m.set_zero();
    CHECK(m.is_zero(0, 0));
    m.add(0, 0, 2.5);
    CHECK_DOUBLE(m.to_double(0, 0), 2.5);
}

// --- large / strided coverage ----------------------------------------------

// A deterministic value per (row, column, batch), spanning a wide exponent
// range with mixed signs so columns rescale and widen at different rates and
// entries partially cancel.
double
wide_sample(std::size_t i, std::size_t j, int batch)
{
    const int e = static_cast<int>((i * 7 + j * 13 + batch * 3) % 120) - 60;
    const double frac =
        static_cast<double>((i * 31 + j * 17 + batch * 11) % 97) / 97.0;
    const double m = 1.0 + frac;
    return std::ldexp((batch % 2 == 0) ? m : -m, e);
}

void
test_row_stride()
{
    const std::size_t rows = 5, cols = 3, stride = cols + 4;

    // The padding is NaN, which add() rejects: if any of it is ever read as
    // matrix data, these tests throw rather than quietly passing.
    std::vector<double> padded(rows * stride,
                               std::numeric_limits<double>::quiet_NaN());
    std::vector<double> packed(rows * cols);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            const double v = wide_sample(i, j, 0);
            padded[i * stride + j] = v;
            packed[i * cols + j] = v;
        }
    }

    cbfp::ColumnBlockMatrix a(rows, cols);
    cbfp::ColumnBlockMatrix b(rows, cols);
    a.reserve_for(padded.data(), 4, stride);
    b.reserve_for(packed.data(), 4);

    for (int n = 0; n < 3; ++n) {
        a.add_matrix(padded.data(), stride);
        b.add_matrix(packed.data());
    }
    a.add_matrix_scaled_pow2(padded.data(), -5, stride);
    b.add_matrix_scaled_pow2(packed.data(), -5);

    for (std::size_t j = 0; j < cols; ++j) {
        CHECK(a.column_exponent(j) == b.column_exponent(j));
        CHECK(a.column_bit_width(j) == b.column_bit_width(j));
    }
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            CHECK_STR(a.to_exact_decimal(i, j), b.to_exact_decimal(i, j));
        }
    }

    // Reading back into a padded buffer must leave the padding untouched.
    const double kFill = -12345.0;
    std::vector<double> out(rows * stride, kFill);
    a.to_matrix(out.data(), stride);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            CHECK_DOUBLE(out[i * stride + j], b.to_double(i, j));
        }
        for (std::size_t j = cols; j < stride; ++j) {
            CHECK_DOUBLE(out[i * stride + j], kFill);
        }
    }
}

void
test_large_matrix_matches_scalar()
{
    // A column of 257 rows shares one exponent across every row, while a 1x1
    // accumulator picks the exponent that suits its single cell. The stored
    // integers therefore differ, but the exact values must not: this pins down
    // both the row-stride arithmetic and the claim that a column's shared
    // scale never changes what it holds.
    const std::size_t rows = 257, cols = 9;
    const int batches = 8;

    cbfp::ColumnBlockMatrix big(rows, cols);
    std::vector<double> buf(rows * cols);
    for (int b = 0; b < batches; ++b) {
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                buf[i * cols + j] = wide_sample(i, j, b);
            }
        }
        big.add_matrix(buf.data());
    }

    // Every column had to rescale and widen well past its first value.
    for (std::size_t j = 0; j < cols; ++j) {
        CHECK(big.column_bit_width(j) >= 128);
    }

    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            cbfp::ColumnBlockMatrix one(1, 1);
            for (int b = 0; b < batches; ++b) {
                one.add(0, 0, wide_sample(i, j, b));
            }
            CHECK_STR(big.to_exact_decimal(i, j), one.to_exact_decimal(0, 0));
            CHECK_DOUBLE(big.to_double(i, j), one.to_double(0, 0));
            // The shared column scale is at least as fine as the cell needs.
            CHECK(big.column_exponent(j) <= one.column_exponent(0));
        }
    }
}

void
test_large_matrix_cancels_to_zero()
{
    // Subtracting the same values back in a different order must clear every
    // one of the 2313 cells, which no amount of cross-row bleed would survive.
    const std::size_t rows = 257, cols = 9;
    const int batches = 6;

    cbfp::ColumnBlockMatrix m(rows, cols);
    std::vector<double> buf(rows * cols);
    for (int b = 0; b < batches; ++b) {
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                buf[i * cols + j] = wide_sample(i, j, b);
            }
        }
        m.add_matrix(buf.data());
    }

    // Walk the batches back in reverse, and the cells within each batch in
    // reverse too, so nothing is undone in the order it was applied.
    for (int b = batches; b-- > 0;) {
        for (std::size_t i = rows; i-- > 0;) {
            for (std::size_t j = cols; j-- > 0;) {
                m.sub(i, j, wide_sample(i, j, b));
            }
        }
    }

    std::size_t nonzero = 0;
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            if (!m.is_zero(i, j)) ++nonzero;
        }
    }
    CHECK(nonzero == 0);
}

void
test_large_matrix_per_column_scales()
{
    // 512 x 64 with a distinct power-of-two scale per column. Each cell sums
    // integers below 2^53 at a single fixed scale, so plain double addition is
    // itself exact here and serves as an independent reference for all 32768
    // cells -- this test is about indexing and scale bookkeeping at size.
    const std::size_t rows = 512, cols = 64;
    const int batches = 16;

    cbfp::ColumnBlockMatrix m(rows, cols);
    std::vector<double> buf(rows * cols);
    std::vector<double> reference(rows * cols, 0.0);

    for (int b = 0; b < batches; ++b) {
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                const int n =
                    static_cast<int>((i * 13 + j * 7 + b * 5) % 2001) - 1000;
                const double v = std::ldexp(static_cast<double>(n),
                                            static_cast<int>(j) - 32);
                buf[i * cols + j] = v;
                reference[i * cols + j] += v;
            }
        }
        m.add_matrix(buf.data());
    }

    std::vector<double> out(rows * cols, 0.0);
    m.to_matrix(out.data());

    std::size_t mismatches = 0;
    for (std::size_t k = 0; k < rows * cols; ++k) {
        if (out[k] != reference[k]) ++mismatches;
    }
    CHECK(mismatches == 0);

    // Integer values at a fixed per-column scale stay narrow: the column
    // exponent should land on the column's own power of two, not drift down.
    for (std::size_t j = 0; j < cols; ++j) {
        CHECK(m.column_exponent(j) >= static_cast<int>(j) - 32);
        CHECK(m.column_bit_width(j) <= 128);
    }
}

void
test_column_major_input()
{
    // The same matrix in both layouts must accumulate to the same exact
    // values. Column-major needs no staging copy, so it takes a different
    // path through the kernels.
    const std::size_t rows = 133, cols = 7;
    const int batches = 5;

    std::vector<double> row_major(rows * cols);
    std::vector<double> col_major(rows * cols);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            const double v = wide_sample(i, j, 0);
            row_major[i * cols + j] = v;
            col_major[j * rows + i] = v;
        }
    }

    cbfp::ColumnBlockMatrix a(rows, cols);
    cbfp::ColumnBlockMatrix b(rows, cols);
    for (int n = 0; n < batches; ++n) {
        a.add_matrix(row_major.data());
        b.add_matrix_col_major(col_major.data());
    }
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            CHECK_STR(a.to_exact_decimal(i, j), b.to_exact_decimal(i, j));
        }
    }
    for (std::size_t j = 0; j < cols; ++j) {
        CHECK(a.column_exponent(j) == b.column_exponent(j));
        CHECK(a.column_bit_width(j) == b.column_bit_width(j));
    }

    // Power-of-two scaling agrees too.
    cbfp::ColumnBlockMatrix c(rows, cols);
    cbfp::ColumnBlockMatrix d(rows, cols);
    c.add_matrix_scaled_pow2(row_major.data(), -7);
    d.add_matrix_col_major_scaled_pow2(col_major.data(), -7);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            CHECK_STR(c.to_exact_decimal(i, j), d.to_exact_decimal(i, j));
        }
    }

    // A padded column stride must leave the gaps unread: NaN there would be
    // rejected if it were ever treated as data.
    const std::size_t pad = rows + 5;
    std::vector<double> padded(pad * cols,
                               std::numeric_limits<double>::quiet_NaN());
    for (std::size_t j = 0; j < cols; ++j) {
        for (std::size_t i = 0; i < rows; ++i) {
            padded[j * pad + i] = col_major[j * rows + i];
        }
    }
    cbfp::ColumnBlockMatrix e(rows, cols);
    for (int n = 0; n < batches; ++n)
        e.add_matrix_col_major(padded.data(), pad);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            CHECK_STR(e.to_exact_decimal(i, j), a.to_exact_decimal(i, j));
        }
    }
}

void
test_streamed_columns()
{
    // Feeding one column at a time must match feeding the whole matrix. This
    // is what lets a caller hold an accumulator far larger than any input it
    // could keep resident.
    const std::size_t rows = 97, cols = 11;
    const int batches = 4;

    std::vector<double> col_major(rows * cols);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            col_major[j * rows + i] = wide_sample(i, j, 0);
        }
    }

    cbfp::ColumnBlockMatrix whole(rows, cols);
    cbfp::ColumnBlockMatrix streamed(rows, cols);
    for (int n = 0; n < batches; ++n) {
        whole.add_matrix_col_major(col_major.data());
        // Deliberately out of order: columns are independent.
        for (std::size_t k = 0; k < cols; ++k) {
            const std::size_t j = (k * 7 + 3) % cols;
            streamed.add_column(j, col_major.data() + j * rows);
        }
    }
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            CHECK_STR(whole.to_exact_decimal(i, j),
                      streamed.to_exact_decimal(i, j));
        }
    }
    for (std::size_t j = 0; j < cols; ++j) {
        CHECK(whole.column_exponent(j) == streamed.column_exponent(j));
        CHECK(whole.column_bit_width(j) == streamed.column_bit_width(j));
    }

    // Scaled form agrees too, and a bad index is rejected.
    cbfp::ColumnBlockMatrix a(rows, cols), b(rows, cols);
    a.add_matrix_col_major_scaled_pow2(col_major.data(), -9);
    for (std::size_t j = 0; j < cols; ++j) {
        b.add_column_scaled_pow2(j, col_major.data() + j * rows, -9);
    }
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            CHECK_STR(a.to_exact_decimal(i, j), b.to_exact_decimal(i, j));
        }
    }

    bool threw = false;
    try {
        b.add_column(cols, col_major.data());
    } catch (const std::out_of_range &) {
        threw = true;
    }
    CHECK(threw);
}

// Value for matrix `m` at (i, j). Each matrix occupies a distinct magnitude
// band, so different orderings trigger the column rescales at different points
// and pass through different intermediate widths on the way to the same total.
double
order_sample(std::size_t i, std::size_t j, int m)
{
    const int band = (m * 137) % 400 - 200;
    const int e = band + static_cast<int>((i * 5 + j * 3) % 17);
    const double frac =
        1.0 + static_cast<double>((i * 31 + j * 17 + m * 7) % 101) / 101.0;
    return std::ldexp((m % 3 == 0) ? -frac : frac, e);
}

void
test_order_independence()
{
    // Double addition is already commutative -- a + b and b + a round
    // identically -- but it is not associative, which is what makes a naive
    // running total depend on arrival order. Exact accumulation restores
    // associativity, and the two together give order independence: reordering
    // a running total needs both, since
    //     (a+b)+c = a+(b+c) = a+(c+b) = (a+c)+b
    // is assoc, comm, assoc. So any permutation must agree digit for digit.
    const std::size_t rows = 23, cols = 7;
    const int count = 10;

    std::vector<std::vector<double>> mats;
    for (int m = 0; m < count; ++m) {
        std::vector<double> a(rows * cols);
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                a[i * cols + j] = order_sample(i, j, m);
            }
        }
        mats.push_back(std::move(a));
    }

    cbfp::ColumnBlockMatrix ref(rows, cols);
    for (int m = 0; m < count; ++m) ref.add_matrix(mats[m].data());

    std::vector<std::string> want(rows * cols);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            want[i * cols + j] = ref.to_exact_decimal(i, j);
        }
    }

    std::mt19937_64 rng(31337);
    std::vector<int> perm(count);
    for (int m = 0; m < count; ++m) perm[m] = m;

    for (int trial = 0; trial < 24; ++trial) {
        std::shuffle(perm.begin(), perm.end(), rng);
        cbfp::ColumnBlockMatrix t(rows, cols);
        for (int k = 0; k < count; ++k) t.add_matrix(mats[perm[k]].data());

        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                CHECK_STR(t.to_exact_decimal(i, j), want[i * cols + j]);
            }
        }
        // The representation lands in the same place too, which is a stronger
        // claim than associativity and does not follow from it: the exponent
        // is a minimum over every value seen and the width a maximum over the
        // same set, neither depending on arrival order. It is not guaranteed
        // in general -- a column that passes through exact zero lets rescale
        // reset the width bound -- but it holds wherever that does not occur.
        for (std::size_t j = 0; j < cols; ++j) {
            CHECK(t.column_exponent(j) == ref.column_exponent(j));
            CHECK(t.column_bit_width(j) == ref.column_bit_width(j));
        }
    }
}

void
test_order_independent_cancellation()
{
    // Matrices that cancel exactly must reach zero in any order, including
    // orders where the huge terms arrive before the tiny ones and orders where
    // the running total passes through values far larger than the answer.
    const std::size_t rows = 17, cols = 5;
    const int pairs = 6;

    std::vector<std::vector<double>> mats;
    for (int m = 0; m < pairs; ++m) {
        std::vector<double> a(rows * cols), b(rows * cols);
        for (std::size_t k = 0; k < rows * cols; ++k) {
            a[k] = order_sample(k / cols, k % cols, m);
            b[k] = -a[k];
        }
        mats.push_back(std::move(a));
        mats.push_back(std::move(b));
    }

    std::mt19937_64 rng(4242);
    std::vector<int> perm(mats.size());
    for (std::size_t m = 0; m < mats.size(); ++m) perm[m] = static_cast<int>(m);

    for (int trial = 0; trial < 20; ++trial) {
        std::shuffle(perm.begin(), perm.end(), rng);
        cbfp::ColumnBlockMatrix t(rows, cols);
        for (int k : perm) t.add_matrix(mats[k].data());

        std::size_t nonzero = 0;
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                if (!t.is_zero(i, j)) ++nonzero;
            }
        }
        CHECK(nonzero == 0);
    }
}

void
test_order_independent_across_entry_points()
{
    // Whole-matrix, column-major, per-column and per-element accumulation are
    // four different code paths; all must agree on the same total.
    const std::size_t rows = 29, cols = 6;
    const int count = 8;

    std::vector<std::vector<double>> row_major, col_major;
    for (int m = 0; m < count; ++m) {
        std::vector<double> r(rows * cols), c(rows * cols);
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                const double v = order_sample(i, j, m);
                r[i * cols + j] = v;
                c[j * rows + i] = v;
            }
        }
        row_major.push_back(std::move(r));
        col_major.push_back(std::move(c));
    }

    cbfp::ColumnBlockMatrix by_matrix(rows, cols);
    cbfp::ColumnBlockMatrix by_col_major(rows, cols);
    cbfp::ColumnBlockMatrix by_column(rows, cols);
    cbfp::ColumnBlockMatrix by_element(rows, cols);

    std::mt19937_64 rng(909);
    std::vector<int> perm(count);
    for (int m = 0; m < count; ++m) perm[m] = m;

    for (int m = 0; m < count; ++m) by_matrix.add_matrix(row_major[m].data());

    std::shuffle(perm.begin(), perm.end(), rng);
    for (int m : perm) by_col_major.add_matrix_col_major(col_major[m].data());

    std::shuffle(perm.begin(), perm.end(), rng);
    for (int m : perm) {
        for (std::size_t j = 0; j < cols; ++j) {
            by_column.add_column(j, col_major[m].data() + j * rows);
        }
    }

    std::shuffle(perm.begin(), perm.end(), rng);
    for (int m : perm) {
        for (std::size_t j = cols; j-- > 0;) {
            for (std::size_t i = rows; i-- > 0;) {
                const double v = order_sample(i, j, m);
                // Exercise sub() as well: adding -v and subtracting v must be
                // the same operation.
                if ((i + j) & 1) {
                    by_element.sub(i, j, -v);
                } else {
                    by_element.add(i, j, v);
                }
            }
        }
    }

    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            const std::string want = by_matrix.to_exact_decimal(i, j);
            CHECK_STR(by_col_major.to_exact_decimal(i, j), want);
            CHECK_STR(by_column.to_exact_decimal(i, j), want);
            CHECK_STR(by_element.to_exact_decimal(i, j), want);
        }
    }
}

void
test_zero_crossing_preserves_value()
{
    // A column that passes through exact zero lets rescale take its shortcut
    // and reset the width bound, so two orderings can end with different
    // widths. The value must be identical regardless: the width bounds what
    // can be stored, it is never part of what is stored.
    const std::size_t rows = 11, cols = 3;

    auto fill = [rows, cols](int base) {
        std::vector<double> a(rows * cols);
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                a[i * cols + j] = std::ldexp(
                    1.0 + static_cast<double>((i * 7 + j * 5) % 37) / 37.0,
                    base + static_cast<int>((i + j) % 5));
            }
        }
        return a;
    };

    const std::vector<double> big = fill(100);
    std::vector<double> neg_big = big;
    for (double &x : neg_big) x = -x;             // cancels `big` exactly
    const std::vector<double> tiny = fill(-200);  // forces a deep rescale
    const std::vector<double> mid = fill(40);

    // Order A: the column reaches exact zero, and only then does a rescale
    // arrive -- so it takes the all-zero shortcut and resets the bound.
    cbfp::ColumnBlockMatrix a(rows, cols);
    a.add_matrix(big.data());
    a.add_matrix(neg_big.data());
    a.add_matrix(tiny.data());
    a.add_matrix(mid.data());

    // Order B: the same values, but the rescale happens while the column
    // still holds data, so the bound is carried across it instead.
    cbfp::ColumnBlockMatrix b(rows, cols);
    b.add_matrix(big.data());
    b.add_matrix(tiny.data());
    b.add_matrix(neg_big.data());
    b.add_matrix(mid.data());

    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            CHECK_STR(a.to_exact_decimal(i, j), b.to_exact_decimal(i, j));
            CHECK_DOUBLE(a.to_double(i, j), b.to_double(i, j));
        }
    }
    for (std::size_t j = 0; j < cols; ++j) {
        CHECK(a.column_exponent(j) == b.column_exponent(j));
    }

    // Whether the widths diverge is incidental; report it so the exception
    // stays visible rather than becoming folklore.
    std::size_t widest_a = 0, widest_b = 0;
    for (std::size_t j = 0; j < cols; ++j) {
        widest_a = std::max(widest_a, a.column_bit_width(j));
        widest_b = std::max(widest_b, b.column_bit_width(j));
    }
    std::printf("  [zero-crossing widths: %zu vs %zu bits, values identical]\n",
                widest_a, widest_b);
}

}  // namespace

// Threading partitions by column, and columns are independent in storage, so a
// threaded run must produce not merely the same value but the same bits: the
// same exponent, the same width, the same limbs. Anything else would mean the
// partitioning had leaked.
static void
test_threaded_matches_serial()
{
    const std::size_t rows = 257, cols = 33;
    std::mt19937_64 rng(90210);
    const int nbatches = 4;
    std::vector<std::vector<double>> batch(nbatches,
                                           std::vector<double>(rows * cols));
    for (auto &b : batch) {
        for (auto &x : b) {
            const std::uint64_t m = (rng() | (std::uint64_t{1} << 52)) &
                                    ((std::uint64_t{1} << 53) - 1);
            x = std::ldexp(static_cast<double>(m),
                           static_cast<int>(rng() % 300) - 150);
            if (rng() & 1) x = -x;
        }
    }

    cbfp::ColumnBlockMatrix serial(rows, cols);
    for (const auto &b : batch) serial.add_matrix_col_major(b.data());

    for (unsigned n : {2u, 4u, 8u}) {
        cbfp::ColumnBlockMatrix threaded(rows, cols);
        threaded.set_threads(n);
        CHECK(threaded.threads() == n);
        for (const auto &b : batch) threaded.add_matrix_col_major(b.data());

        std::size_t differing = 0;
        for (std::size_t j = 0; j < cols; ++j) {
            if (threaded.column_exponent(j) != serial.column_exponent(j) ||
                threaded.column_limbs(j) != serial.column_limbs(j)) {
                ++differing;
                continue;
            }
            for (std::size_t i = 0; i < rows; ++i) {
                if (threaded.entry_limbs(i, j) != serial.entry_limbs(i, j)) {
                    ++differing;
                }
            }
        }
        CHECK(0 == differing);  // bit-identical to the serial run
    }

    // The row-major path stages each column through a buffer; with threads
    // that buffer has to be per worker, which is the one piece of state the
    // columns used to share.
    {
        std::vector<double> rowmajor(rows * cols);
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                rowmajor[i * cols + j] = batch[0][j * rows + i];
            }
        }
        cbfp::ColumnBlockMatrix one(rows, cols), many(rows, cols);
        many.set_threads(8);
        one.add_matrix(rowmajor.data());
        many.add_matrix(rowmajor.data());
        std::size_t differing = 0;
        for (std::size_t j = 0; j < cols; ++j) {
            for (std::size_t i = 0; i < rows; ++i) {
                if (one.entry_limbs(i, j) != many.entry_limbs(i, j))
                    ++differing;
            }
        }
        CHECK(0 == differing);  // per-worker staging, not shared
    }
}

// The residual is the exact difference between the stored value and the double
// to_double returns. Checking that without a second arbitrary-precision
// implementation means making the container check itself: accumulate the same
// values again, subtract the rounded result, and whatever remains must be
// exactly what the residual reported.
static void
test_readback_residual()
{
    const std::size_t rows = 48, cols = 4;
    const int reps = 8;
    std::mt19937_64 rng(90210);

    // Distinct batches, not one batch added repeatedly: adding the same value
    // eight times only multiplies it by eight, which is exact and would leave
    // every residual zero.
    std::vector<std::vector<double>> batches(reps,
                                             std::vector<double>(rows * cols));
    for (auto &batch : batches) {
        for (auto &x : batch) {
            const std::uint64_t m = (rng() | (std::uint64_t{1} << 52)) &
                                    ((std::uint64_t{1} << 53) - 1);
            // A wide exponent spread, so a cell's terms do not share a scale
            // and the sum needs far more than 53 bits to hold exactly.
            x = std::ldexp(static_cast<double>(m),
                           static_cast<int>(rng() % 80) - 40);
            if (rng() & 1) x = -x;
        }
    }

    cbfp::ColumnBlockMatrix a(rows, cols);
    for (int r = 0; r < reps; ++r) a.add_matrix(batches[r].data());

    int nonzero_residuals = 0, negative_residuals = 0;
    for (std::size_t j = 0; j < cols; ++j) {
        for (std::size_t i = 0; i < rows; ++i) {
            double lo = 1.0;  // must be overwritten, never left as-is
            const double hi = a.to_double(i, j, &lo);

            // Asking for the residual must not perturb the value.
            CHECK_DOUBLE(a.to_double(i, j), hi);

            // The exact statement of what a residual is: re-accumulate the
            // same value, subtract the double that was returned, and what is
            // left must round to the residual.
            cbfp::ColumnBlockMatrix b(1, 1);
            for (int r = 0; r < reps; ++r) {
                b.add_column(0, &batches[r][i * cols + j]);
            }
            const double minus_hi = -hi;
            b.add_column(0, &minus_hi);
            CHECK_DOUBLE(b.to_double(0, 0), lo);

            // A two-term expansion: the residual never reaches a full ulp of
            // the value it corrects, so the two do not overlap.
            if (0.0 != hi) {
                const double ulp =
                    std::nextafter(std::fabs(hi),
                                   std::numeric_limits<double>::infinity()) -
                    std::fabs(hi);
                CHECK(std::fabs(lo) <= ulp / 2);
            }

            // Exactly representable and "nothing was discarded" must be the
            // same predicate.
            CHECK((0.0 == lo) == a.is_exactly_representable(i, j));

            if (0.0 != lo) ++nonzero_residuals;
            if (lo < 0.0) ++negative_residuals;
        }
    }
    // The test is only meaningful if it actually exercised rounding in both
    // directions rather than landing on exact values throughout.
    CHECK(nonzero_residuals > static_cast<int>(rows * cols / 2));
    CHECK(negative_residuals > 0);

    // An exactly representable entry reports a residual of zero, positive.
    cbfp::ColumnBlockMatrix e(1, 1);
    const double two = 2.0;
    e.add_column(0, &two);
    double elo = 1.0;
    CHECK_DOUBLE(e.to_double(0, 0, &elo), 2.0);
    CHECK_DOUBLE(elo, 0.0);
    CHECK(!std::signbit(elo));

    // So does an empty one.
    cbfp::ColumnBlockMatrix z(1, 1);
    double zlo = 1.0;
    CHECK_DOUBLE(z.to_double(0, 0, &zlo), 0.0);
    CHECK_DOUBLE(zlo, 0.0);

    // The bulk form must agree with the per-entry form, stride included.
    const std::size_t stride = cols + 3;
    std::vector<double> mv(rows * stride, -1.0), mr(rows * stride, -1.0);
    a.to_matrix_with_residual(mv.data(), mr.data(), stride);
    for (std::size_t j = 0; j < cols; ++j) {
        for (std::size_t i = 0; i < rows; ++i) {
            double want_lo = 0.0;
            const double want_hi = a.to_double(i, j, &want_lo);
            CHECK_DOUBLE(mv[i * stride + j], want_hi);
            CHECK_DOUBLE(mr[i * stride + j], want_lo);
        }
    }
}

// The classic case the container exists for: a sum whose exact value needs far
// more than a double, where the residual says how much the answer was off by.
static void
test_residual_recovers_cancellation()
{
    cbfp::ColumnBlockMatrix a(1, 1);
    const double big = 1e300, one = 1.0;
    a.add_column(0, &big);
    a.add_column(0, &one);

    // 1e300 + 1 is not representable, so the double is just 1e300 and the
    // residual is what naive summation silently threw away.
    double lo = 0.0;
    const double hi = a.to_double(0, 0, &lo);
    CHECK_DOUBLE(hi, 1e300);
    CHECK_DOUBLE(lo, 1.0);

    // 0.1 added ten times. The double 0.1 is 3602879701896397 * 2^-55, so ten
    // of them come to 4503599627370496.25 * 2^-52 -- a quarter of the way
    // above 1.0, which rounds down to exactly 1.0 and leaves 2^-54 behind.
    cbfp::ColumnBlockMatrix t(1, 1);
    const double tenth = 0.1;
    for (int k = 0; k < 10; ++k) t.add_column(0, &tenth);
    double tlo = 0.0;
    const double thi = t.to_double(0, 0, &tlo);
    CHECK_DOUBLE(thi, 1.0);
    CHECK_DOUBLE(tlo, std::ldexp(1.0, -54));

    // And the pair reconstructs: subtracting the double leaves the residual.
    cbfp::ColumnBlockMatrix u(1, 1);
    for (int k = 0; k < 10; ++k) u.add_column(0, &tenth);
    const double minus = -thi;
    u.add_column(0, &minus);
    CHECK_DOUBLE(u.to_double(0, 0), tlo);
}

// Storing one triangle has to be invisible in the answers. The reference is a
// full-storage accumulator fed the same symmetric matrices: every entry must
// agree, on both sides of the diagonal and for either uplo.
static void
test_symmetric_matches_full_storage()
{
    const std::size_t n = 37;
    const int reps = 5;
    std::mt19937_64 rng(5150);

    std::vector<std::vector<double>> batches(reps, std::vector<double>(n * n));
    for (auto &b : batches) {
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j <= i; ++j) {
                const std::uint64_t m = (rng() | (std::uint64_t{1} << 52)) &
                                        ((std::uint64_t{1} << 53) - 1);
                double v = std::ldexp(static_cast<double>(m),
                                      static_cast<int>(rng() % 60) - 30);
                if (rng() & 1) v = -v;
                b[i * n + j] = v;
                b[j * n + i] = v;  // symmetric by construction
            }
        }
    }

    cbfp::ColumnBlockMatrix full(n, n);
    cbfp::ColumnBlockMatrix lo(n, cbfp::Uplo::Lower);
    cbfp::ColumnBlockMatrix up(n, cbfp::Uplo::Upper);
    for (int r = 0; r < reps; ++r) {
        full.add_matrix(batches[r].data());
        lo.add_matrix(batches[r].data());
        up.add_matrix(batches[r].data());
    }

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            const double want = full.to_double(i, j);
            CHECK_DOUBLE(lo.to_double(i, j), want);
            CHECK_DOUBLE(up.to_double(i, j), want);
        }
    }

    // Stored once, so the two halves are the same limbs and not merely the
    // same value -- which is what full storage cannot promise, since columns
    // i and j carry their own exponents.
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            CHECK(lo.entry_limbs(i, j) == lo.entry_limbs(j, i));
            CHECK(up.entry_limbs(i, j) == up.entry_limbs(j, i));
        }
    }

    // The shape itself.
    CHECK(lo.symmetric() && !full.symmetric());
    for (std::size_t j = 0; j < n; ++j) {
        CHECK(lo.column_rows(j) == n - j);
        CHECK(lo.column_first_row(j) == j);
        CHECK(up.column_rows(j) == j + 1);
        CHECK(up.column_first_row(j) == 0);
    }

    // The limb payload halves, but memory_bytes() also charges a fixed cost
    // per limb column, and at small n that cost dominates: the measured ratio
    // is 0.85 at n=37, 0.71 at 128, 0.57 at 512 and 0.52 at 2048. So check the
    // saving where the payload actually dominates rather than asserting a
    // bound the small case cannot meet.
    CHECK(lo.memory_bytes() < full.memory_bytes());
    {
        const std::size_t big = 512;
        cbfp::ColumnBlockMatrix bf(big, big);
        cbfp::ColumnBlockMatrix bl(big, cbfp::Uplo::Lower);
        cbfp::ColumnBlockMatrix bu(big, cbfp::Uplo::Upper);
        CHECK(bl.memory_bytes() < bf.memory_bytes() * 3 / 4);
        // The two triangles are the same multiset of column lengths, so they
        // cost exactly the same however the padding falls.
        CHECK(bl.memory_bytes() == bu.memory_bytes());
    }

    // to_matrix fills the whole n x n, mirrored.
    std::vector<double> out(n * n, -1.0);
    lo.to_matrix(out.data());
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            CHECK_DOUBLE(out[i * n + j], full.to_double(i, j));
            CHECK_DOUBLE(out[i * n + j], out[j * n + i]);
        }
    }
}

// The column-major entry point, the threaded path and the weighted partition
// all have to reach the same answer as the serial row-major one.
static void
test_symmetric_entry_points_agree()
{
    const std::size_t n = 53;
    std::mt19937_64 rng(2718);

    std::vector<double> rowmajor(n * n), colmajor(n * n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            const std::uint64_t m = (rng() | (std::uint64_t{1} << 52)) &
                                    ((std::uint64_t{1} << 53) - 1);
            double v = std::ldexp(static_cast<double>(m),
                                  static_cast<int>(rng() % 40) - 20);
            if (rng() & 1) v = -v;
            rowmajor[i * n + j] = rowmajor[j * n + i] = v;
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) colmajor[j * n + i] = rowmajor[i * n + j];
    }

    cbfp::ColumnBlockMatrix a(n, cbfp::Uplo::Lower);
    a.add_matrix(rowmajor.data());

    cbfp::ColumnBlockMatrix b(n, cbfp::Uplo::Lower);
    b.add_matrix_col_major(colmajor.data());

    cbfp::ColumnBlockMatrix c(n, cbfp::Uplo::Lower);
    c.set_threads(4);
    c.add_matrix(rowmajor.data());

    // Streamed one column at a time: each takes only its stored slice.
    cbfp::ColumnBlockMatrix d(n, cbfp::Uplo::Lower);
    for (std::size_t j = 0; j < n; ++j) {
        d.add_column(j, &colmajor[j * n + d.column_first_row(j)]);
    }

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            const double want = a.to_double(i, j);
            CHECK_DOUBLE(b.to_double(i, j), want);
            CHECK_DOUBLE(c.to_double(i, j), want);
            CHECK_DOUBLE(d.to_double(i, j), want);
            CHECK(a.entry_limbs(i, j) == c.entry_limbs(i, j));
        }
    }

    // The upper triangle stored column-major is a different slice, so check it
    // separately rather than assuming the lower case covered it.
    cbfp::ColumnBlockMatrix e(n, cbfp::Uplo::Upper);
    e.add_matrix_col_major(colmajor.data());
    cbfp::ColumnBlockMatrix f(n, cbfp::Uplo::Upper);
    f.set_threads(3);
    f.add_matrix(rowmajor.data());
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            CHECK_DOUBLE(e.to_double(i, j), a.to_double(i, j));
            CHECK_DOUBLE(f.to_double(i, j), a.to_double(i, j));
        }
    }
}

// A symmetric accumulator still has to do everything the general one does:
// widen, rescale, report exactly, and pre-size from surveys.
static void
test_symmetric_keeps_exactness()
{
    // The diagonal is stored once and must be added once, not twice.
    cbfp::ColumnBlockMatrix d(3, cbfp::Uplo::Lower);
    std::vector<double> m(9, 0.0);
    m[0] = 2.0;   // (0,0)
    m[4] = -7.5;  // (1,1)
    d.add_matrix(m.data());
    d.add_matrix(m.data());
    CHECK_DOUBLE(d.to_double(0, 0), 4.0);
    CHECK_DOUBLE(d.to_double(1, 1), -15.0);
    CHECK_DOUBLE(d.to_double(2, 2), 0.0);

    // Cancellation across the whole double range, on an off-diagonal entry.
    cbfp::ColumnBlockMatrix a(4, cbfp::Uplo::Lower);
    a.add(2, 1, 1e300);
    a.add(2, 1, 1.0);
    a.add(1, 2, -1e300);  // the mirrored index reaches the same entry
    CHECK_DOUBLE(a.to_double(2, 1), 1.0);
    CHECK_DOUBLE(a.to_double(1, 2), 1.0);
    CHECK_STR(a.to_exact_decimal(1, 2), "1");

    // sub through the mirror cancels exactly.
    a.sub(1, 2, 1.0);
    CHECK(a.is_zero(2, 1));

    // A value far below the column's scale forces a rescale of a triangular
    // column; nothing is lost and the mirror still agrees.
    cbfp::ColumnBlockMatrix r(6, cbfp::Uplo::Upper);
    r.add(4, 2, 1.0);
    r.add(4, 2, std::ldexp(1.0, -300));
    CHECK(r.column_exponent(4) == -300);
    r.sub(2, 4, 1.0);
    CHECK_DOUBLE(r.to_double(2, 4), std::ldexp(1.0, -300));
    CHECK(r.is_exactly_representable(4, 2));

    // Residuals work through the mirror too.
    cbfp::ColumnBlockMatrix q(3, cbfp::Uplo::Lower);
    for (int k = 0; k < 10; ++k) q.add(2, 0, 0.1);
    double lo = 0.0;
    CHECK_DOUBLE(q.to_double(0, 2, &lo), 1.0);
    CHECK_DOUBLE(lo, std::ldexp(1.0, -54));

    // Pre-sized from surveys of the stored triangle: same answers, and no
    // rescale or widen during accumulation.
    const std::size_t n = 16;
    std::mt19937_64 rng(99);
    std::vector<double> b(n * n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double v = std::ldexp(1.0 + static_cast<double>(rng() % 1000),
                                  static_cast<int>(rng() % 20) - 10);
            b[i * n + j] = b[j * n + i] = v;
        }
    }
    cbfp::ColumnBlockMatrix adaptive(n, cbfp::Uplo::Lower);
    adaptive.add_matrix(b.data());
    adaptive.add_matrix(b.data());

    std::vector<std::vector<cbfp::Survey>> surveys(
        2, std::vector<cbfp::Survey>(n));
    for (std::size_t j = 0; j < n; ++j) {
        std::vector<double> slice;
        for (std::size_t k = 0; k < n - j; ++k) slice.push_back(b[(j + k) * n + j]);
        surveys[0][j] = surveys[1][j] =
            cbfp::survey_column(slice.data(), slice.size());
    }
    cbfp::ColumnBlockMatrix presized(n, cbfp::Uplo::Lower, surveys);
    const std::size_t limbs_before = presized.column_limbs(0);
    const int exp_before = presized.column_exponent(0);
    presized.add_matrix(b.data());
    presized.add_matrix(b.data());
    CHECK(presized.column_limbs(0) == limbs_before);
    CHECK(presized.column_exponent(0) == exp_before);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            CHECK_DOUBLE(presized.to_double(i, j), adaptive.to_double(i, j));
        }
    }
}

// Folding several matrices into one pass changes only the traffic. Exact
// accumulation does not care about order or grouping, so the answer has to be
// identical to adding them one at a time.
static void
test_fold_matches_sequential()
{
    const std::size_t rows = 137, cols = 7;
    std::mt19937_64 rng(60613);

    // Three spreads: narrow enough for one limb, wide enough for several, and
    // wide enough to exceed the fold's register budget and fall back.
    for (int spread : {6, 90, 900}) {
        const int nbatches = 5;
        std::vector<std::vector<double>> batches(
            nbatches, std::vector<double>(rows * cols));
        for (auto &b : batches) {
            for (auto &x : b) {
                const std::uint64_t m = (rng() | (std::uint64_t{1} << 52)) &
                                        ((std::uint64_t{1} << 53) - 1);
                x = std::ldexp(static_cast<double>(m),
                               static_cast<int>(rng() % (spread + 1)) -
                                   spread / 2);
                if (rng() & 1) x = -x;
            }
        }
        std::vector<const double *> ptrs;
        for (auto &b : batches) ptrs.push_back(b.data());

        cbfp::ColumnBlockMatrix seq(rows, cols);
        for (auto &b : batches) seq.add_matrix_col_major(b.data());

        cbfp::ColumnBlockMatrix fold(rows, cols);
        fold.add_matrices_col_major(ptrs.data(), ptrs.size());

        // Exact decimal rather than limbs: an adaptive container that rescaled
        // partway through can end up allocated wider than one told the extents
        // up front, which changes the limbs without changing the value.
        std::size_t bad = 0;
        for (std::size_t j = 0; j < cols; ++j) {
            for (std::size_t i = 0; i < rows; ++i) {
                if (seq.to_exact_decimal(i, j) != fold.to_exact_decimal(i, j)) {
                    ++bad;
                }
            }
        }
        CHECK(0 == bad);
        // The widest spread must actually have exercised the fallback.
        if (900 == spread) CHECK(fold.column_limbs(0) > 8);
        if (6 == spread) CHECK(fold.column_limbs(0) <= 8);
    }

    // Pre-sized on both sides, where the reservation is identical, so the
    // stored limbs must match bit for bit and not merely the value.
    {
        const int nbatches = 6;
        std::vector<std::vector<double>> batches(
            nbatches, std::vector<double>(rows * cols));
        for (auto &b : batches) {
            for (auto &x : b) {
                const std::uint64_t m = (rng() | (std::uint64_t{1} << 52)) &
                                        ((std::uint64_t{1} << 53) - 1);
                x = std::ldexp(static_cast<double>(m),
                               static_cast<int>(rng() % 21) - 10);
                if (rng() & 1) x = -x;
            }
        }
        std::vector<const double *> ptrs;
        for (auto &b : batches) ptrs.push_back(b.data());

        std::vector<std::vector<cbfp::Survey>> surveys;
        for (auto &b : batches) {
            std::vector<cbfp::Survey> s(cols);
            for (std::size_t j = 0; j < cols; ++j) {
                s[j] = cbfp::survey_column(b.data() + j * rows, rows);
            }
            surveys.push_back(s);
        }

        cbfp::ColumnBlockMatrix seq(rows, cols, surveys);
        for (auto &b : batches) seq.add_matrix_col_major(b.data());

        cbfp::ColumnBlockMatrix fold(rows, cols, surveys);
        fold.add_matrices_col_major(ptrs.data(), ptrs.size());

        std::size_t bad = 0;
        for (std::size_t j = 0; j < cols; ++j) {
            if (seq.column_exponent(j) != fold.column_exponent(j) ||
                seq.column_limbs(j) != fold.column_limbs(j)) {
                ++bad;
                continue;
            }
            for (std::size_t i = 0; i < rows; ++i) {
                if (seq.entry_limbs(i, j) != fold.entry_limbs(i, j)) ++bad;
            }
        }
        CHECK(0 == bad);

        // And the declared count is spent by the fold, not by the call.
        std::vector<const double *> one{ptrs[0]};
        bool threw = false;
        try {
            fold.add_matrices_col_major(one.data(), 1);
        } catch (const std::runtime_error &) {
            threw = true;
        }
        CHECK(threw);
    }

    // Threaded, symmetric, and a degenerate count of one.
    {
        const std::size_t n = 64;
        std::vector<std::vector<double>> batches(3, std::vector<double>(n * n));
        for (auto &b : batches) {
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t j = 0; j <= i; ++j) {
                    const std::uint64_t m = (rng() | (std::uint64_t{1} << 52)) &
                                            ((std::uint64_t{1} << 53) - 1);
                    double v = std::ldexp(static_cast<double>(m),
                                          static_cast<int>(rng() % 13) - 6);
                    if (rng() & 1) v = -v;
                    b[j * n + i] = b[i * n + j] = v;
                }
            }
        }
        std::vector<const double *> ptrs;
        for (auto &b : batches) ptrs.push_back(b.data());

        cbfp::ColumnBlockMatrix seq(n, cbfp::Uplo::Lower);
        for (auto &b : batches) seq.add_matrix_col_major(b.data());

        cbfp::ColumnBlockMatrix fold(n, cbfp::Uplo::Lower);
        fold.set_threads(4);
        fold.add_matrices_col_major(ptrs.data(), ptrs.size());

        std::size_t bad = 0;
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                if (seq.to_exact_decimal(i, j) != fold.to_exact_decimal(i, j)) {
                    ++bad;
                }
            }
        }
        CHECK(0 == bad);

        cbfp::ColumnBlockMatrix one(n, cbfp::Uplo::Lower);
        one.add_matrices_col_major(ptrs.data(), 1);
        cbfp::ColumnBlockMatrix plain(n, cbfp::Uplo::Lower);
        plain.add_matrix_col_major(batches[0].data());
        bad = 0;
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                if (one.entry_limbs(i, j) != plain.entry_limbs(i, j)) ++bad;
            }
        }
        CHECK(0 == bad);
    }

    // A fold that contradicts its metadata still reports.
    {
        std::vector<double> v(64, 0.5);  // exponent -1
        const double *p[1] = {v.data()};
        std::vector<std::vector<cbfp::Survey>> lie{
            {cbfp::Survey{0, 8, true, false}}};
        cbfp::ColumnBlockMatrix m(64, 1, lie);
        bool threw = false;
        try {
            m.add_matrices_col_major(p, 1);
        } catch (const std::runtime_error &) {
            threw = true;
        }
        CHECK(threw);
    }
}

// The public survey is what a producer computes so an accumulator need not.
// It has to agree with what the accumulator would have worked out itself,
// or preallocating from it is worse than useless.
static void
test_public_survey()
{
    const std::size_t rows = 133, cols = 5;
    std::mt19937_64 rng(31337);
    std::vector<double> colmajor(rows * cols), rowmajor(rows * cols);
    for (std::size_t j = 0; j < cols; ++j) {
        for (std::size_t i = 0; i < rows; ++i) {
            const std::uint64_t m = (rng() | (std::uint64_t{1} << 52)) &
                                    ((std::uint64_t{1} << 53) - 1);
            double v = std::ldexp(static_cast<double>(m),
                                  static_cast<int>(rng() % 200) - 100);
            if (rng() & 1) v = -v;
            if (0 == j) v = 0.0;  // an all-zero column
            colmajor[j * rows + i] = v;
            rowmajor[i * cols + j] = v;
        }
    }

    std::vector<cbfp::Survey> a(cols), b(cols);
    cbfp::survey_matrix_col_major(colmajor.data(), rows, cols, a.data());
    cbfp::survey_matrix(rowmajor.data(), rows, cols, b.data());

    for (std::size_t j = 0; j < cols; ++j) {
        CHECK(a[j].any == b[j].any);
        CHECK(a[j].nonfinite == b[j].nonfinite);
        CHECK(a[j].min_exponent == b[j].min_exponent);
        CHECK(a[j].max_top == b[j].max_top);
    }
    CHECK(!a[0].any);  // the zero column
    CHECK(a[1].any);

    // The extents must actually bound the data: accumulating one matrix into a
    // column reserved from its own survey must need no rescale and no widen.
    for (std::size_t j = 1; j < cols; ++j) {
        cbfp::ColumnBlockMatrix one(rows, 1);
        one.reserve_column(
            j == 1 ? 0 : 0, static_cast<int>(a[j].min_exponent),
            static_cast<std::size_t>(a[j].max_top - a[j].min_exponent) + 2);
        const int exp_before = one.column_exponent(0);
        const std::size_t limbs_before = one.column_limbs(0);
        one.add_column(0, colmajor.data() + j * rows);
        CHECK(one.column_exponent(0) == exp_before);
        CHECK(one.column_limbs(0) == limbs_before);
    }

    // A non-finite anywhere must be reported.
    {
        std::vector<double> bad(rows, 1.0);
        bad[rows / 2] = std::numeric_limits<double>::infinity();
        const cbfp::Survey s = cbfp::survey_column(bad.data(), rows);
        CHECK(s.nonfinite);
    }
}

// An accumulator given the surveys of every matrix it will receive must reach
// exactly the same limbs as one that worked them out for itself -- otherwise
// pre-sizing is a different algorithm rather than the same one told in advance.
static void
test_presized_matches_adaptive()
{
    const std::size_t rows = 211, cols = 9;
    const int nbatches = 6;
    std::mt19937_64 rng(56789);
    std::vector<std::vector<double>> batch(nbatches,
                                           std::vector<double>(rows * cols));
    for (auto &b : batch) {
        for (auto &x : b) {
            const std::uint64_t m = (rng() | (std::uint64_t{1} << 52)) &
                                    ((std::uint64_t{1} << 53) - 1);
            x = std::ldexp(static_cast<double>(m),
                           static_cast<int>(rng() % 260) - 130);
            if (rng() & 1) x = -x;
        }
    }

    // What a producer would hand forward.
    std::vector<std::vector<cbfp::Survey>> surveys(
        nbatches, std::vector<cbfp::Survey>(cols));
    for (int b = 0; b < nbatches; ++b) {
        cbfp::survey_matrix_col_major(batch[b].data(), rows, cols,
                                      surveys[b].data());
    }

    cbfp::ColumnBlockMatrix adaptive(rows, cols);
    cbfp::ColumnBlockMatrix presized(rows, cols, surveys);
    for (int b = 0; b < nbatches; ++b) {
        adaptive.add_matrix_col_major(batch[b].data());
        presized.add_matrix_col_major(batch[b].data());
    }

    for (std::size_t j = 0; j < cols; ++j) {
        CHECK(presized.column_exponent(j) == adaptive.column_exponent(j));
        for (std::size_t i = 0; i < rows; ++i) {
            // Widths may differ -- the adaptive one grows as it learns -- so
            // compare the value, sign-extending the narrower.
            const std::vector<std::uint64_t> a = adaptive.entry_limbs(i, j);
            const std::vector<std::uint64_t> p = presized.entry_limbs(i, j);
            const std::size_t n = std::max(a.size(), p.size());
            const std::uint64_t fa =
                0 != (a.back() >> 63) ? ~std::uint64_t{0} : 0;
            const std::uint64_t fp =
                0 != (p.back() >> 63) ? ~std::uint64_t{0} : 0;
            for (std::size_t k = 0; k < n; ++k) {
                CHECK((k < a.size() ? a[k] : fa) == (k < p.size() ? p[k] : fp));
            }
        }
    }

    // Submitting in a different order must land in the same place, since the
    // reservation is an aggregate and knows nothing about which matrix is
    // which.
    {
        cbfp::ColumnBlockMatrix reversed(rows, cols, surveys);
        for (int b = nbatches; b-- > 0;) {
            reversed.add_matrix_col_major(batch[b].data());
        }
        for (std::size_t j = 0; j < cols; ++j) {
            for (std::size_t i = 0; i < rows; ++i) {
                CHECK(reversed.entry_limbs(i, j) == presized.entry_limbs(i, j));
            }
        }
    }

    // A matrix outside what was declared is detected and reported rather than
    // quietly producing a wrong sum. The accumulator is not usable afterwards,
    // which is the point: the contract was broken, not accommodated.
    {
        std::vector<std::vector<cbfp::Survey>> one(
            1, std::vector<cbfp::Survey>(1));
        std::vector<double> modest(rows, 1.0);
        cbfp::survey_matrix_col_major(modest.data(), rows, 1, one[0].data());

        // ...then hand it something far outside those extents.
        std::vector<double> huge(rows, std::ldexp(1.0, 400));
        cbfp::ColumnBlockMatrix lied(rows, 1, one);
        bool threw = false;
        try {
            lied.add_matrix_col_major(huge.data());
        } catch (const std::runtime_error &) {
            threw = true;
        }
        CHECK(threw);

        // And a value below the declared exponent, which is the other
        // direction and the one that silently dropped values before.
        std::vector<double> tiny(rows, std::ldexp(1.0, -400));
        cbfp::ColumnBlockMatrix lied2(rows, 1, one);
        threw = false;
        try {
            lied2.add_matrix_col_major(tiny.data());
        } catch (const std::runtime_error &) {
            threw = true;
        }
        CHECK(threw);
    }

    // The count is the one part of the contract that is enforced, because it
    // costs nothing: the width bound holds for that many matrices, not more.
    {
        cbfp::ColumnBlockMatrix full(rows, cols, surveys);
        for (int b = 0; b < nbatches; ++b)
            full.add_matrix_col_major(batch[0].data());
        bool threw = false;
        try {
            full.add_matrix_col_major(batch[0].data());
        } catch (const std::runtime_error &) {
            threw = true;
        }
        CHECK(threw);
    }
}

int
main()
{
    test_decompose();
    test_single_value_roundtrip();
    test_catastrophic_cancellation();
    test_repeated_tenth();
    test_exact_decimal();
    test_exponent_and_width_tracking();
    test_full_double_range();
    test_rounding_ties_to_even();
    test_add_then_subtract_is_zero();
    test_matrix_accumulation();
    test_reserve_avoids_rescaling();
    test_column_independence();
    test_errors();
    test_many_small_into_large();
    test_reuse_after_zero();
    test_row_stride();
    test_large_matrix_matches_scalar();
    test_large_matrix_cancels_to_zero();
    test_large_matrix_per_column_scales();
    test_column_major_input();
    test_streamed_columns();
    test_order_independence();
    test_order_independent_cancellation();
    test_order_independent_across_entry_points();
    test_zero_crossing_preserves_value();
    test_threaded_matches_serial();
    test_public_survey();
    test_presized_matches_adaptive();
    test_readback_residual();
    test_residual_recovers_cancellation();
    test_symmetric_matches_full_storage();
    test_symmetric_entry_points_agree();
    test_symmetric_keeps_exactness();
    test_fold_matches_sequential();

    std::printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
