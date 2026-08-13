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
#include <cuda_runtime.h>

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
check(bool ok, const std::string &what)
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
cross_check(const char *name, std::size_t rows, std::size_t cols, int nbatches,
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
    for (auto &v : first) v = random_value(200);
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
            for (auto &one : b) {
                for (auto &x : one) x = random_value(spread);
            }
            cbfp::ColumnBlockMatrix cpu(rows, cols);
            for (auto &one : b) cpu.add_matrix_col_major(one.data());
            cbfp::CudaColumnBlockMatrix gpu(rows, cols);
            gpu.reserve_like(cpu);
            for (auto &one : b) gpu.add_matrix_col_major(one.data());

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

// reserve_for is what lets the device container stand on its own. Its contract
// is the CPU's: pre-size from a representative batch about to be accumulated
// `count` times, and accumulation then never needs to grow a column.
void
test_reserve_for()
{
    const std::size_t rows = 251, cols = 7;
    const int nbatches = 5;
    std::vector<double> batch(rows * cols);
    for (auto &x : batch) x = random_value(120);
    // one column entirely zero: it must still be reserved, since an
    // unreserved column is an error on the device rather than something that
    // can grow on first use
    for (std::size_t i = 0; i < rows; ++i) batch[3 * rows + i] = 0.0;

    cbfp::ColumnBlockMatrix cpu(rows, cols);
    cbfp::CudaColumnBlockMatrix gpu(rows, cols);
    gpu.reserve_for(batch.data(), nbatches);

    bool threw = false;
    try {
        for (int b = 0; b < nbatches; ++b) {
            cpu.add_matrix_col_major(batch.data());
            gpu.add_matrix_col_major(batch.data());
        }
    } catch (const std::exception &) {
        threw = true;
    }
    check(!threw, "reserve_for holds for the count it was given");

    std::size_t mismatches = 0;
    for (std::size_t j = 0; j < cols; ++j) {
        const std::vector<std::uint64_t> dev = gpu.download_column(j);
        for (std::size_t i = 0; i < rows; ++i) {
            const std::vector<std::uint64_t> host = cpu.entry_limbs(i, j);
            // The two may hold different widths -- the CPU grows as it goes,
            // reserve_for sizes up front -- so compare the value, not the
            // representation, over the limbs they share plus sign extension.
            const std::size_t n = gpu.column_limbs(j);
            if (gpu.column_exponent(j) != cpu.column_exponent(j)) {
                ++mismatches;
                break;
            }
            for (std::size_t k = 0; k < n; ++k) {
                const std::uint64_t h =
                    k < host.size()
                        ? host[k]
                        : (0 != (host.back() >> 63) ? ~std::uint64_t{0} : 0);
                if (h != dev[k * rows + i]) ++mismatches;
            }
        }
    }
    check(0 == mismatches,
          "reserve_for gives bit-identical values to the CPU (" +
              std::to_string(mismatches) + " limb mismatches)");

    // Surveying the same batch on the device must reserve identically.
    {
        double *dev = nullptr;
        cudaMalloc(&dev, batch.size() * sizeof(double));
        cudaMemcpy(dev, batch.data(), batch.size() * sizeof(double),
                   cudaMemcpyHostToDevice);
        cbfp::CudaColumnBlockMatrix g2(rows, cols);
        g2.reserve_for_device(dev, nbatches);
        std::size_t differ = 0;
        for (std::size_t j = 0; j < cols; ++j) {
            if (g2.column_exponent(j) != gpu.column_exponent(j) ||
                g2.column_limbs(j) != gpu.column_limbs(j)) {
                ++differ;
            }
        }
        check(0 == differ,
              "reserve_for_device reserves the same as reserve_for");
        cudaFree(dev);
    }
}

// An addend does not grow with the batch count but the accumulated sum does,
// so checking only the incoming values against the reservation is not enough:
// twenty batches of 2^60 into a one-limb column once wrapped modulo 2^64 and
// reported nothing. The column now grows to fit instead, which in a limb-major
// layout is an allocation and a sign fill with no data movement.
void
test_accumulated_bound()
{
    const std::size_t rows = 8, cols = 1;
    const int nbatches = 20;
    std::vector<double> b(rows * cols);
    for (std::size_t i = 0; i < rows; ++i) b[i] = std::ldexp(1.0, 60);

    cbfp::ColumnBlockMatrix cpu(rows, cols);
    cbfp::CudaColumnBlockMatrix gpu(rows, cols);
    gpu.reserve_column(0, 0, 64);  // one limb, deliberately too narrow

    for (int k = 0; k < nbatches; ++k) {
        gpu.add_matrix_col_major(b.data());
        cpu.add_matrix_col_major(b.data());
    }
    check(gpu.column_limbs(0) > 1, "the column grew past its reservation");

    // The CPU chose exponent 60 and the device was pinned to 0, so the device
    // holds the same value scaled by 2^60. Compare that, not the encoding.
    std::size_t mismatches = 0;
    for (std::size_t i = 0; i < rows; ++i) {
        const std::vector<std::uint64_t> host = cpu.entry_limbs(i, 0);
        const std::vector<std::uint64_t> dev = gpu.entry_limbs(i, 0);
        unsigned __int128 want = 0;
        for (std::size_t k = host.size(); k-- > 0;) {
            want = (want << 64) | host[k];
        }
        want <<= 60;
        for (std::size_t k = 0; k < dev.size(); ++k) {
            const std::uint64_t w =
                static_cast<std::uint64_t>(want >> (64 * k));
            if (w != dev[k]) ++mismatches;
        }
    }
    check(0 == mismatches, "growing kept the sum exact (" +
                               std::to_string(mismatches) +
                               " limb mismatches)");
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
        } catch (const std::domain_error &) {
            threw = true;
        }
        check(threw, "NaN is rejected by the device survey");

        bad[70] = std::numeric_limits<double>::infinity();
        threw = false;
        try {
            gpu.add_matrix_col_major(bad.data());
        } catch (const std::domain_error &) {
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
        } catch (const std::runtime_error &) {
            threw = true;
        }
        check(threw, "an under-reserved exponent is an error");
    }

    // A width that is too small is no longer an error -- the column grows.
    // Only the exponent cannot be fixed after the fact, since that needs every
    // entry shifted rather than a limb appended.
    {
        cbfp::CudaColumnBlockMatrix gpu(64, 1);
        gpu.reserve_column(0, 0, 64);
        std::vector<double> v(64, std::ldexp(1.0, 100));
        bool threw = false;
        try {
            gpu.add_matrix_col_major(v.data());
        } catch (const std::runtime_error &) {
            threw = true;
        }
        check(!threw, "an under-reserved width grows instead of failing");
        check(gpu.column_limbs(0) > 1, "and the column is wider for it");
    }

    // Accumulating into a column that was never reserved.
    {
        cbfp::CudaColumnBlockMatrix gpu(64, 2);
        gpu.reserve_column(0, 0, 64);
        std::vector<double> v(64 * 2, 1.0);
        bool threw = false;
        try {
            gpu.add_matrix_col_major(v.data());
        } catch (const std::runtime_error &) {
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
    test_reserve_for();
    test_accumulated_bound();
    test_errors();

    std::printf("%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
