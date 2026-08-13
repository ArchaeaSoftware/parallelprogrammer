// Cross-checks the device accumulator against the CPU one, bit for bit.
//
// The two hold the same integer only if they agree about every step: the
// IEEE-754 field extraction, the odd normalization, the limb offset, the
// two-limb addend split, and the carry chain. Comparing stored limbs is
// strictly stronger than comparing to_double or to_exact_decimal, since
// identical limbs at an identical column exponent imply identical everything
// downstream -- and unlike a decimal comparison it cannot be passed by two
// implementations that are wrong in the same way at readback.
//
// Same assertion-runner style as tests/test_cbfp.cpp, no dependencies.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "cbfp/column_accumulator.hpp"
#include "cbfp/cuda_accumulator.hpp"

namespace {

int checks = 0;
int failures = 0;

void
check(bool ok, const std::string& what)
{
    ++checks;
    if (!ok) {
        ++failures;
        std::printf("FAIL: %s\n", what.c_str());
    }
}

// Runs the same column-major batches through both implementations and
// compares every stored limb.
//
// The CPU container learns each column's exponent and width while
// accumulating; the device one is then reserved to exactly that, which is the
// only way the two can be compared. Rescaling is exact, so a column that
// reached its final exponent by rescaling holds the same integer as one that
// started there.
template <typename Gen>
void
cross_check(const char* name, std::size_t rows, std::size_t cols, int nbatches,
            Gen gen)
{
    std::vector<std::vector<double>> batches(
        nbatches, std::vector<double>(rows * cols, 0.0));
    for (int b = 0; b < nbatches; ++b) {
        for (std::size_t j = 0; j < cols; ++j) {
            for (std::size_t i = 0; i < rows; ++i) {
                batches[b][j * rows + i] = gen(i, j, b);
            }
        }
    }

    cbfp::ColumnBlockMatrix cpu(rows, cols);
    for (int b = 0; b < nbatches; ++b) {
        cpu.add_matrix_col_major(batches[b].data());
    }

    cbfp::CudaColumnBlockMatrix gpu(rows, cols);
    gpu.reserve_like(cpu);
    for (int b = 0; b < nbatches; ++b) {
        gpu.add_matrix_col_major(batches[b].data());
    }

    std::size_t mismatches = 0;
    std::string first;
    for (std::size_t j = 0; j < cols; ++j) {
        if (gpu.column_exponent(j) != cpu.column_exponent(j) ||
            gpu.column_limbs(j) != cpu.column_limbs(j)) {
            ++mismatches;
            continue;
        }
        const std::vector<std::uint64_t> dev = gpu.download_column(j);
        const std::size_t n = cpu.column_limbs(j);
        for (std::size_t i = 0; i < rows; ++i) {
            const std::vector<std::uint64_t> host = cpu.entry_limbs(i, j);
            for (std::size_t k = 0; k < n; ++k) {
                if (host[k] == dev[k * rows + i]) continue;
                ++mismatches;
                if (first.empty()) {
                    char buf[256];
                    std::snprintf(buf, sizeof buf,
                                  " first at (%zu,%zu) limb %zu: cpu %016llx "
                                  "gpu %016llx",
                                  i, j, k, (unsigned long long)host[k],
                                  (unsigned long long)dev[k * rows + i]);
                    first = buf;
                }
            }
        }
    }
    check(mismatches == 0, std::string(name) + ": " +
                               std::to_string(mismatches) + " mismatches" +
                               first);
}

std::mt19937_64 rng(20260812);

// A random finite double whose exponent lands within `spread` binades.
double
random_value(int spread)
{
    std::uint64_t m =
        (rng() | (std::uint64_t{1} << 52)) & ((std::uint64_t{1} << 53) - 1);
    const int e = spread == 0 ? 0 : static_cast<int>(rng() % (spread + 1));
    const double v = std::ldexp(static_cast<double>(m), e - spread / 2);
    return (rng() & 1) ? v : -v;
}

void
test_shapes()
{
    // The workhorse shape, exponents clustered.
    cross_check("4096x64 narrow, 4 batches", 4096, 64, 4,
                [](std::size_t, std::size_t, int) { return random_value(8); });

    // A ragged row count, which the grid-stride loop has to tail correctly.
    // 257 is the shape the CPU suite uses for the same reason.
    cross_check("257x9 ragged, 3 batches", 257, 9, 3,
                [](std::size_t, std::size_t, int) { return random_value(16); });

    // Another ragged count, one row short of a warp multiple.
    cross_check("1021x5 ragged, 2 batches", 1021, 5, 2,
                [](std::size_t, std::size_t, int) { return random_value(32); });

    // Divergent exponents: threads in a warp land in different limb arrays,
    // and each walks a different limb range.
    cross_check(
        "512x8 wide spread, 3 batches", 512, 8, 3,
        [](std::size_t, std::size_t, int) { return random_value(400); });

    // A single column and a single row, where the grid degenerates.
    cross_check("4096x1 single column", 4096, 1, 2,
                [](std::size_t, std::size_t, int) { return random_value(24); });
    cross_check("1x64 single row", 1, 64, 2,
                [](std::size_t, std::size_t, int) { return random_value(24); });
}

void
test_edge_values()
{
    // Zeros, negative zero, subnormals and the extremes of the range, mixed
    // with ordinary values so columns still have a sensible scale.
    const double specials[] = {0.0,
                               -0.0,
                               1.0,
                               -1.0,
                               0x1p-1074,  // smallest subnormal
                               -0x1p-1074,
                               0x1.fffffffffffffp-1023,  // largest subnormal
                               0x1p-1022,                // smallest normal
                               0x1.fffffffffffffp+52};
    const std::size_t n = sizeof(specials) / sizeof(specials[0]);

    cross_check("specials cycled, 3 batches", 64, 9, 3,
                [&](std::size_t i, std::size_t j, int b) {
                    return specials[(i + j + static_cast<std::size_t>(b)) % n];
                });

    // Cancellation: the second batch is the negation of the first, so every
    // cell must return to exact zero on both sides.
    std::vector<double> first(128 * 4);
    for (auto& v : first) v = random_value(200);
    cross_check("cancellation to exact zero", 128, 4, 2,
                [&](std::size_t i, std::size_t j, int b) {
                    const double v = first[j * 128 + i];
                    return b == 0 ? v : -v;
                });
}

// The accumulate reports, per column, the highest limb position holding
// anything but sign extension. The fit test that will consume it needs one
// property above all: nothing significant may live above the reported figure,
// or a column would be judged to fit when it does not.
void
test_occupancy()
{
    // A borrow out of a negative addend writes all-ones every limb to the top
    // of a column. Counting written limbs rather than significant ones would
    // report the full width for any column that has ever gone negative --
    // which is exactly the cancellation-heavy case the figure exists for.
    {
        const std::size_t rows = 64, cols = 3;
        std::vector<double> v(rows * cols, 0.0);
        for (std::size_t i = 0; i < rows; ++i) {
            v[0 * rows + i] = 1.0;
            v[1 * rows + i] = (7 == i) ? -1.0 : 1.0;
            v[2 * rows + i] = -1.0;
        }
        cbfp::CudaColumnBlockMatrix gpu(rows, cols);
        for (std::size_t j = 0; j < cols; ++j) gpu.reserve_column(j, 0, 256);
        gpu.add_matrix_col_major(v.data());
        check(0 == gpu.column_occupancy(0), "occupancy: all positive, limb 0");
        check(0 == gpu.column_occupancy(1),
              "occupancy: one negative does not saturate to nlimbs-1");
        check(-1 == gpu.column_occupancy(2),
              "occupancy: -1 in every limb is pure sign, so nothing is "
              "significant");
    }

    // Nothing significant above the reported occupancy, over random shapes,
    // spreads and batch counts -- including columns wide enough for a
    // high-water mark set by one batch to be missed by a later one.
    {
        std::size_t checked = 0, violations = 0;
        for (int trial = 0; trial < 12; ++trial) {
            const std::size_t rows = 97, cols = 5;
            const int spread = static_cast<int>(rng() % 400);
            const int nbatches = 1 + static_cast<int>(rng() % 4);
            std::vector<std::vector<double>> b(
                nbatches, std::vector<double>(rows * cols));
            for (auto& one : b) {
                for (auto& x : one) x = random_value(spread);
            }
            cbfp::ColumnBlockMatrix cpu(rows, cols);
            for (auto& one : b) cpu.add_matrix_col_major(one.data());
            cbfp::CudaColumnBlockMatrix gpu(rows, cols);
            gpu.reserve_like(cpu);
            for (auto& one : b) gpu.add_matrix_col_major(one.data());

            for (std::size_t j = 0; j < cols; ++j) {
                const int occ = gpu.column_occupancy(j);
                const std::size_t n = cpu.column_limbs(j);
                for (std::size_t i = 0; i < rows; ++i) {
                    const std::vector<std::uint64_t> l = cpu.entry_limbs(i, j);
                    const std::uint64_t fill =
                        0 != (l[n - 1] >> 63) ? ~std::uint64_t{0} : 0;
                    int top = -1;
                    for (int k = static_cast<int>(n) - 1; k >= 0; --k) {
                        if (l[k] != fill) {
                            top = k;
                            break;
                        }
                    }
                    ++checked;
                    if (top > occ) ++violations;
                }
            }
        }
        check(0 == violations, "occupancy bounds every significant limb (" +
                                   std::to_string(checked) + " entries, " +
                                   std::to_string(violations) + " above)");
    }
}

void
test_errors()
{
    // A non-finite value must be rejected by the survey, not silently added.
    {
        cbfp::ColumnBlockMatrix cpu(64, 2);
        std::vector<double> ok(64 * 2, 1.0);
        cpu.add_matrix_col_major(ok.data());

        cbfp::CudaColumnBlockMatrix gpu(64, 2);
        gpu.reserve_like(cpu);
        std::vector<double> bad = ok;
        bad[70] = std::nan("");
        bool threw = false;
        try {
            gpu.add_matrix_col_major(bad.data());
        } catch (const std::domain_error&) {
            threw = true;
        }
        check(threw, "NaN is rejected by the device survey");

        bad[70] = std::numeric_limits<double>::infinity();
        threw = false;
        try {
            gpu.add_matrix_col_major(bad.data());
        } catch (const std::domain_error&) {
            threw = true;
        }
        check(threw, "inf is rejected by the device survey");
    }

    // A column reserved at too high an exponent must be an error rather than
    // a silent loss of precision: rescaling is impossible mid-launch.
    {
        cbfp::CudaColumnBlockMatrix gpu(64, 1);
        gpu.reserve_column(0, 0, 128);  // exponent 2^0, so no fractions fit
        std::vector<double> v(64, 0.5);
        bool threw = false;
        try {
            gpu.add_matrix_col_major(v.data());
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw, "an under-reserved exponent is an error");
    }

    // Likewise a column reserved too narrow to hold what arrives.
    {
        cbfp::CudaColumnBlockMatrix gpu(64, 1);
        gpu.reserve_column(0, 0, 64);
        std::vector<double> v(64, std::ldexp(1.0, 100));
        bool threw = false;
        try {
            gpu.add_matrix_col_major(v.data());
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw, "an under-reserved width is an error");
    }

    // Accumulating into a column that was never reserved.
    {
        cbfp::CudaColumnBlockMatrix gpu(64, 2);
        gpu.reserve_column(0, 0, 64);
        std::vector<double> v(64 * 2, 1.0);
        bool threw = false;
        try {
            gpu.add_matrix_col_major(v.data());
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw, "an unreserved column is an error");
    }
}

}  // namespace

int
main()
{
    if (!cbfp::cuda_available()) {
        std::printf("no CUDA device; skipping\n");
        return 0;
    }

    test_shapes();
    test_edge_values();
    test_occupancy();
    test_errors();

    std::printf("%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
