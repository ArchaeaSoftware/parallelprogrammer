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
#include <cstring>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
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

// The CPU container still takes ordinary memory, so this is a pass-through.
// Having both overloads lets a cross-check read the same on either side.
void
submit(cbfp::ColumnBlockMatrix &m, const std::vector<double> &v,
       std::size_t col_stride = 0)
{
    m.add_matrix_col_major(v.data(), col_stride);
}

// The device path takes page-locked input only, so this is the copy a caller
// now has to write for themselves -- made once, here, where it is visible.
// It waits on the handle before the pinned buffer goes out of scope, which is
// exactly the contract the API states rather than the library guessing at it.
void
submit(cbfp::CudaColumnBlockMatrix &g, const std::vector<double> &v,
       std::size_t col_stride = 0)
{
    double *p = nullptr;
    if (cudaSuccess !=
        cudaHostAlloc(&p, v.size() * sizeof(double), cudaHostAllocMapped)) {
        throw std::runtime_error("test: cudaHostAlloc failed");
    }
    std::memcpy(p, v.data(), v.size() * sizeof(double));
    try {
        cbfp::InputRead h = g.add_matrix_col_major(p, col_stride);
        h.wait();
    } catch (...) {
        cudaFreeHost(p);
        throw;
    }
    cudaFreeHost(p);
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
        submit(cpu, batches[b]);
    }

    cbfp::CudaColumnBlockMatrix gpu(rows, cols);
    gpu.reserve_like(cpu);
    for (int b = 0; b < nbatches; ++b) {
        submit(gpu, batches[b]);
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
        submit(gpu, v);
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
            for (auto &one : b) submit(cpu, one);
            cbfp::CudaColumnBlockMatrix gpu(rows, cols);
            gpu.reserve_like(cpu);
            for (auto &one : b) submit(gpu, one);

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
            submit(cpu, batch);
            submit(gpu, batch);
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
        submit(gpu, b);
        submit(cpu, b);
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

// Compares one entry by value across two containers that may have settled on
// different exponents and widths: the device holds V * 2^egpu, the CPU
// V * 2^ecpu, so the CPU's limbs shifted left by the difference must match.
bool
same_value(const cbfp::ColumnBlockMatrix &cpu,
           const cbfp::CudaColumnBlockMatrix &gpu, std::size_t i, std::size_t j)
{
    const std::vector<std::uint64_t> h = cpu.entry_limbs(i, j);
    const std::vector<std::uint64_t> d = gpu.entry_limbs(i, j);
    const long ecpu = cpu.column_exponent(j);
    const long egpu = gpu.column_exponent(j);
    if (egpu > ecpu) return false;  // the device must never end up coarser
    const unsigned long sh = static_cast<unsigned long>(ecpu - egpu);
    const bool neg = 0 != (h.back() >> 63);

    const unsigned long word = sh / 64;
    const unsigned bit = static_cast<unsigned>(sh % 64);
    for (std::size_t k = 0; k < d.size(); ++k) {
        const long src = static_cast<long>(k) - static_cast<long>(word);
        auto at = [&](long q) -> std::uint64_t {
            if (q < 0) return 0;
            if (q < static_cast<long>(h.size())) return h[q];
            return neg ? ~std::uint64_t{0} : 0;
        };
        std::uint64_t want;
        if (src < 0) {
            want = 0;
        } else if (0 == bit) {
            want = at(src);
        } else {
            want = (at(src) << bit) | (at(src - 1) >> (64 - bit));
        }
        if (want != d[k]) return false;
    }
    return true;
}

// Lowering a column's exponent shifts every entry rather than appending to
// them, so it is the half that cannot be repaired inside a launch. Done
// between launches it must still land on exactly the value the CPU holds.
void
test_rescale()
{
    const std::size_t rows = 131, cols = 4;
    std::mt19937_64 r(4242);

    // Batches whose ulps descend, so each one forces the exponent lower.
    const int nbatches = 5;
    std::vector<std::vector<double>> b(nbatches,
                                       std::vector<double>(rows * cols));
    for (int k = 0; k < nbatches; ++k) {
        for (auto &x : b[k]) {
            std::uint64_t m = (r() | (std::uint64_t{1} << 52)) &
                              ((std::uint64_t{1} << 53) - 1);
            x = std::ldexp(static_cast<double>(m), 40 - 20 * k);
            if (r() & 1) x = -x;
        }
    }

    cbfp::ColumnBlockMatrix cpu(rows, cols);
    cbfp::CudaColumnBlockMatrix gpu(rows, cols);
    // Deliberately sized for the first batch only, so every later batch
    // forces both a rescale and a widen.
    gpu.reserve_for(b[0].data(), 1);

    for (int k = 0; k < nbatches; ++k) {
        submit(cpu, b[k]);
        submit(gpu, b[k]);
    }

    std::size_t coarse = 0, wrong = 0;
    for (std::size_t j = 0; j < cols; ++j) {
        if (gpu.column_exponent(j) > cpu.column_exponent(j)) ++coarse;
        for (std::size_t i = 0; i < rows; ++i) {
            if (!same_value(cpu, gpu, i, j)) ++wrong;
        }
    }
    check(0 == coarse,
          "rescale reached an exponent at least as low as the CPU");
    check(0 == wrong, "rescale kept every entry exact (" +
                          std::to_string(wrong) + " of " +
                          std::to_string(rows * cols) + " wrong)");

    // A shift that is a whole number of limbs, and one that is not, both
    // exercise the in-place walk differently.
    for (unsigned sh : {64u, 3u, 130u}) {
        cbfp::ColumnBlockMatrix c2(8, 1);
        cbfp::CudaColumnBlockMatrix g2(8, 1);
        std::vector<double> hi(8), lo(8);
        for (std::size_t i = 0; i < 8; ++i) {
            hi[i] = std::ldexp(1.0 + i, 0);
            lo[i] = std::ldexp(1.0 + i, -static_cast<int>(sh));
        }
        submit(c2, hi);
        g2.reserve_for(hi.data(), 2);
        submit(g2, hi);
        submit(c2, lo);
        submit(g2, lo);
        std::size_t bad = 0;
        for (std::size_t i = 0; i < 8; ++i) {
            if (!same_value(c2, g2, i, 0)) ++bad;
        }
        check(0 == bad, "rescale by " + std::to_string(sh) + " bits is exact");
    }
}

// Buffers borrowed from acquire_input are read in place by the kernel rather
// than copied, so the accumulator has to know when the device has finished
// with one before handing it back. Refilling a buffer the device is still
// streaming corrupted about 63% of entries when the fast path was keyed on
// the pointer's memory type instead of on ownership.
void
test_acquired_input()
{
    const std::size_t rows = 4096, cols = 16;
    const int nbatches = 6;
    const std::size_t n = rows * cols;

    std::vector<double> data(n);
    for (auto &x : data) x = random_value(150);

    cbfp::ColumnBlockMatrix cpu(rows, cols);
    for (int b = 0; b < nbatches; ++b) submit(cpu, data);

    // Staged: the caller's buffer is copied, so it may be reused at once.
    cbfp::CudaColumnBlockMatrix staged(rows, cols);
    staged.reserve_like(cpu);
    {
        std::vector<double> buf(n);
        for (int b = 0; b < nbatches; ++b) {
            buf = data;
            submit(staged, buf);
            std::fill(buf.begin(), buf.end(), 0.0);  // safe: it was copied
        }
    }

    // Borrowed: read in place, and acquire_input blocks until that is over.
    // The loop refills as fast as it can, which is the shape that raced.
    cbfp::CudaColumnBlockMatrix borrowed(rows, cols);
    borrowed.reserve_like(cpu);
    for (int b = 0; b < nbatches; ++b) {
        double *buf = borrowed.acquire_input();
        for (std::size_t i = 0; i < n; ++i) buf[i] = data[i];
        borrowed.add_matrix_col_major(buf);
    }

    std::size_t staged_wrong = 0, borrowed_wrong = 0;
    for (std::size_t j = 0; j < cols; ++j) {
        for (std::size_t i = 0; i < rows; ++i) {
            const std::vector<std::uint64_t> h = cpu.entry_limbs(i, j);
            if (staged.entry_limbs(i, j) != h) ++staged_wrong;
            if (borrowed.entry_limbs(i, j) != h) ++borrowed_wrong;
        }
    }
    check(0 == staged_wrong,
          "a staged buffer may be reused as soon as the call returns");
    check(0 == borrowed_wrong,
          "an acquired buffer is not handed back while the device reads it");
}

void
test_errors()
{
    // A non-finite value must be rejected by the survey, not silently added.
    {
        cbfp::ColumnBlockMatrix cpu(64, 2);
        std::vector<double> ok(64 * 2, 1.0);
        submit(cpu, ok);

        cbfp::CudaColumnBlockMatrix gpu(64, 2);
        gpu.reserve_like(cpu);
        std::vector<double> bad = ok;
        bad[70] = std::nan("");
        bool threw = false;
        try {
            submit(gpu, bad);
        } catch (const std::domain_error &) {
            threw = true;
        }
        check(threw, "NaN is rejected by the device survey");

        bad[70] = std::numeric_limits<double>::infinity();
        threw = false;
        try {
            submit(gpu, bad);
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
            submit(gpu, v);
        } catch (const std::runtime_error &) {
            threw = true;
        }
        check(!threw, "an under-reserved exponent rescales instead of failing");
        check(gpu.column_exponent(0) < 0,
              "and the column exponent came down to fit");
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
            submit(gpu, v);
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
            submit(gpu, v);
        } catch (const std::runtime_error &) {
            threw = true;
        }
        check(threw, "an unreserved column is an error");
    }
}

// The same cross-check, on a triangular accumulator. Column j stores
// column_rows(j) entries beginning at column_first_row(j), so both the
// download and the host lookup have to be indexed through that rather than by
// the matrix dimension.
template <typename Gen>
void
cross_check_symmetric(const char *name, std::size_t n, cbfp::Uplo uplo,
                      int nbatches, Gen gen)
{
    std::vector<std::vector<double>> batches(nbatches,
                                             std::vector<double>(n * n, 0.0));
    for (int b = 0; b < nbatches; ++b) {
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j <= i; ++j) {
                const double v = gen(i, j, b);
                batches[b][j * n + i] = v;  // column-major
                batches[b][i * n + j] = v;
            }
        }
    }

    cbfp::ColumnBlockMatrix cpu(n, uplo);
    for (int b = 0; b < nbatches; ++b) {
        submit(cpu, batches[b]);
    }

    cbfp::CudaColumnBlockMatrix gpu(n, uplo);
    gpu.reserve_like(cpu);
    for (int b = 0; b < nbatches; ++b) {
        submit(gpu, batches[b]);
    }

    std::size_t mismatches = 0;
    std::string first;
    for (std::size_t j = 0; j < n; ++j) {
        if (gpu.column_exponent(j) != cpu.column_exponent(j) ||
            gpu.column_limbs(j) != cpu.column_limbs(j) ||
            gpu.column_rows(j) != cpu.column_rows(j) ||
            gpu.column_first_row(j) != cpu.column_first_row(j)) {
            ++mismatches;
            continue;
        }
        const std::vector<std::uint64_t> dev = gpu.download_column(j);
        const std::size_t nl = cpu.column_limbs(j);
        const std::size_t cr = cpu.column_rows(j);
        const std::size_t f0 = cpu.column_first_row(j);
        for (std::size_t s = 0; s < cr; ++s) {
            const std::vector<std::uint64_t> host = cpu.entry_limbs(f0 + s, j);
            for (std::size_t k = 0; k < nl; ++k) {
                if (host[k] == dev[k * cr + s]) continue;
                ++mismatches;
                if (first.empty()) {
                    char buf[256];
                    std::snprintf(buf, sizeof buf,
                                  " first at (%zu,%zu) limb %zu: cpu %016llx "
                                  "gpu %016llx",
                                  f0 + s, j, k, (unsigned long long)host[k],
                                  (unsigned long long)dev[k * cr + s]);
                    first = buf;
                }
            }
        }
    }
    check(mismatches == 0, std::string(name) + ": " +
                               std::to_string(mismatches) + " mismatches" +
                               first);
}

void
test_symmetric()
{
    // Both triangles, and a ragged n so the grid-stride tail is exercised on
    // columns whose lengths are not multiples of a warp.
    cross_check_symmetric(
        "sym 512 lower, 3 batches", 512, cbfp::Uplo::Lower, 3,
        [](std::size_t, std::size_t, int) { return random_value(8); });
    cross_check_symmetric(
        "sym 512 upper, 3 batches", 512, cbfp::Uplo::Upper, 3,
        [](std::size_t, std::size_t, int) { return random_value(8); });
    cross_check_symmetric(
        "sym 257 lower ragged, 2 batches", 257, cbfp::Uplo::Lower, 2,
        [](std::size_t, std::size_t, int) { return random_value(16); });
    cross_check_symmetric(
        "sym 129 upper ragged, 2 batches", 129, cbfp::Uplo::Upper, 2,
        [](std::size_t, std::size_t, int) { return random_value(16); });

    // Divergent exponents, which force rescales and widens on columns whose
    // lengths differ -- the case where a per-column row count that was left
    // as the matrix dimension would read or write past the short columns.
    cross_check_symmetric(
        "sym 128 lower wide spread, 3 batches", 128, cbfp::Uplo::Lower, 3,
        [](std::size_t, std::size_t, int) { return random_value(400); });

    // Degenerate sizes, where the last column holds a single entry.
    cross_check_symmetric(
        "sym 1", 1, cbfp::Uplo::Lower, 2,
        [](std::size_t, std::size_t, int) { return random_value(4); });
    cross_check_symmetric(
        "sym 8 upper", 8, cbfp::Uplo::Upper, 2,
        [](std::size_t, std::size_t, int) { return random_value(4); });

    // The shape itself, and half the device memory.
    {
        cbfp::CudaColumnBlockMatrix lo(256, cbfp::Uplo::Lower);
        cbfp::CudaColumnBlockMatrix up(256, cbfp::Uplo::Upper);
        cbfp::CudaColumnBlockMatrix full(256, 256);
        bool shape_ok = lo.symmetric() && !full.symmetric();
        for (std::size_t j = 0; j < 256; ++j) {
            shape_ok = shape_ok && lo.column_rows(j) == 256 - j &&
                       lo.column_first_row(j) == j &&
                       up.column_rows(j) == j + 1 &&
                       up.column_first_row(j) == 0;
        }
        check(shape_ok, "device triangular columns have the right extents");

        cbfp::ColumnBlockMatrix cpu(256, cbfp::Uplo::Lower);
        std::vector<double> v(256 * 256, 1.0);
        submit(cpu, v);
        lo.reserve_like(cpu);
        full.reserve_like(cbfp::ColumnBlockMatrix(256, 256));
        check(lo.memory_bytes() < full.memory_bytes(),
              "and cost less device memory than full storage");
    }

    // reserve_like has to refuse a mismatched shape: a triangular column is a
    // different length, so full-storage widths would size the wrong entries.
    {
        cbfp::ColumnBlockMatrix cpu_full(64, 64);
        std::vector<double> v(64 * 64, 1.0);
        submit(cpu_full, v);
        cbfp::CudaColumnBlockMatrix gpu_sym(64, cbfp::Uplo::Lower);
        bool threw = false;
        try {
            gpu_sym.reserve_like(cpu_full);
        } catch (const std::runtime_error &) {
            threw = true;
        }
        check(threw, "reserve_like rejects a mismatched symmetry");

        cbfp::ColumnBlockMatrix cpu_lo(64, cbfp::Uplo::Lower);
        submit(cpu_lo, v);
        cbfp::CudaColumnBlockMatrix gpu_up(64, cbfp::Uplo::Upper);
        threw = false;
        try {
            gpu_up.reserve_like(cpu_lo);
        } catch (const std::runtime_error &) {
            threw = true;
        }
        check(threw, "reserve_like rejects a mismatched uplo");
    }
}

// Surveys of a column-major batch, one entry per column, taken over exactly
// the slice the accumulator will read.
std::vector<cbfp::Survey>
survey_batch(const std::vector<double> &b, const cbfp::ColumnBlockMatrix &shape)
{
    std::vector<cbfp::Survey> out(shape.cols());
    for (std::size_t j = 0; j < shape.cols(); ++j) {
        out[j] = cbfp::survey_column(
            b.data() + j * shape.rows() + shape.column_first_row(j),
            shape.column_rows(j));
    }
    return out;
}

void
test_presized()
{
    // Pre-sized on both sides from the same surveys. Both containers run the
    // same reservation arithmetic, so they must agree on exponent and width
    // and therefore on every stored limb -- which a comparison against an
    // adaptively sized container could not check, since headroom would leave
    // it at a different width.
    for (int variant = 0; variant < 3; ++variant) {
        const std::size_t rows = variant == 2 ? 257 : 1024;
        const std::size_t cols = variant == 2 ? 9 : 16;
        const int nbatches = 3;

        std::vector<std::vector<double>> batches(
            nbatches, std::vector<double>(rows * cols, 0.0));
        for (int b = 0; b < nbatches; ++b) {
            for (auto &x : batches[b]) x = random_value(variant == 1 ? 300 : 8);
        }

        cbfp::ColumnBlockMatrix shape(rows, cols);
        std::vector<std::vector<cbfp::Survey>> surveys;
        for (int b = 0; b < nbatches; ++b) {
            surveys.push_back(survey_batch(batches[b], shape));
        }

        cbfp::ColumnBlockMatrix cpu(rows, cols, surveys);
        cbfp::CudaColumnBlockMatrix gpu(rows, cols, surveys);
        for (int b = 0; b < nbatches; ++b) {
            submit(cpu, batches[b]);
            submit(gpu, batches[b]);
        }

        std::size_t mismatches = 0;
        for (std::size_t j = 0; j < cols; ++j) {
            if (gpu.column_exponent(j) != cpu.column_exponent(j) ||
                gpu.column_limbs(j) != cpu.column_limbs(j)) {
                ++mismatches;
                continue;
            }
            const std::vector<std::uint64_t> dev = gpu.download_column(j);
            const std::size_t nl = cpu.column_limbs(j);
            for (std::size_t i = 0; i < rows; ++i) {
                const std::vector<std::uint64_t> host = cpu.entry_limbs(i, j);
                for (std::size_t k = 0; k < nl; ++k) {
                    if (host[k] != dev[k * rows + i]) ++mismatches;
                }
            }
        }
        check(mismatches == 0, "presized device matches presized host, variant " +
                                   std::to_string(variant) + ": " +
                                   std::to_string(mismatches) + " mismatches");
    }

    // Symmetric and pre-sized together.
    {
        const std::size_t n = 129;
        std::vector<double> b(n * n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j <= i; ++j) {
                const double v = random_value(12);
                b[j * n + i] = b[i * n + j] = v;
            }
        }
        cbfp::ColumnBlockMatrix shape(n, cbfp::Uplo::Lower);
        std::vector<std::vector<cbfp::Survey>> surveys{survey_batch(b, shape)};

        cbfp::ColumnBlockMatrix cpu(n, cbfp::Uplo::Lower, surveys);
        cbfp::CudaColumnBlockMatrix gpu(n, cbfp::Uplo::Lower, surveys);
        submit(cpu, b);
        submit(gpu, b);

        std::size_t mismatches = 0;
        for (std::size_t j = 0; j < n; ++j) {
            if (gpu.column_exponent(j) != cpu.column_exponent(j) ||
                gpu.column_limbs(j) != cpu.column_limbs(j)) {
                ++mismatches;
                continue;
            }
            const std::vector<std::uint64_t> dev = gpu.download_column(j);
            const std::size_t cr = cpu.column_rows(j);
            const std::size_t f0 = cpu.column_first_row(j);
            for (std::size_t s2 = 0; s2 < cr; ++s2) {
                const std::vector<std::uint64_t> host =
                    cpu.entry_limbs(f0 + s2, j);
                for (std::size_t k = 0; k < cpu.column_limbs(j); ++k) {
                    if (host[k] != dev[k * cr + s2]) ++mismatches;
                }
            }
        }
        check(mismatches == 0, "presized symmetric device matches host: " +
                                   std::to_string(mismatches) + " mismatches");
    }

    // The count is the one thing checked, because it costs nothing.
    {
        std::vector<double> v(64 * 2, 1.0);
        cbfp::ColumnBlockMatrix shape(64, 2);
        std::vector<std::vector<cbfp::Survey>> surveys{survey_batch(v, shape),
                                                       survey_batch(v, shape)};
        cbfp::CudaColumnBlockMatrix gpu(64, 2, surveys);
        submit(gpu, v);
        submit(gpu, v);
        bool threw = false;
        try {
            submit(gpu, v);
        } catch (const std::runtime_error &) {
            threw = true;
        }
        check(threw, "a third matrix past two declared is an error");
    }
}

// Pre-sizing removes the survey, so nothing prevents a batch that contradicts
// the metadata. The accumulate kernel has to notice instead, and the failure
// it is guarding is silent: a negative shift makes the limb offset enormous
// and the carry loop simply never runs, dropping the value without a word.
void
test_device_detector()
{
    const std::size_t rows = 128;

    // An exponent below what the column was sized for.
    {
        std::vector<double> v(rows, 0.5);  // exponent -1
        cbfp::Survey lie{0, 8, true, false};  // claims nothing below 2^0
        std::vector<std::vector<cbfp::Survey>> surveys{{lie}};
        cbfp::CudaColumnBlockMatrix gpu(rows, 1, surveys);
        submit(gpu, v);
        const unsigned f = gpu.column_contradictions(0);
        check(0 != (f & 2u), "a too-low exponent is detected");

        bool threw = false;
        try {
            gpu.synchronize();
        } catch (const std::runtime_error &) {
            threw = true;
        }
        check(threw, "and reported at the next synchronization");
    }

    // An addend past the column's width.
    {
        std::vector<double> v(rows, std::ldexp(1.0, 300));
        cbfp::Survey lie{0, 8, true, false};
        std::vector<std::vector<cbfp::Survey>> surveys{{lie}};
        cbfp::CudaColumnBlockMatrix gpu(rows, 1, surveys);
        submit(gpu, v);
        check(0 != (gpu.column_contradictions(0) & 4u),
              "an addend past the width is detected");
    }

    // A non-finite value, which the survey would have rejected outright.
    {
        std::vector<double> v(rows, 1.0);
        v[rows / 2] = std::numeric_limits<double>::infinity();
        cbfp::Survey ok{0, 8, true, false};
        std::vector<std::vector<cbfp::Survey>> surveys{{ok}};
        cbfp::CudaColumnBlockMatrix gpu(rows, 1, surveys);
        submit(gpu, v);
        check(0 != (gpu.column_contradictions(0) & 1u),
              "a non-finite value is detected");
    }

    // And an honest batch reports nothing, so the detector is not simply
    // always firing.
    {
        std::vector<double> v(rows * 3);
        for (auto &x : v) x = random_value(6);
        cbfp::ColumnBlockMatrix shape(rows, 3);
        std::vector<std::vector<cbfp::Survey>> surveys{survey_batch(v, shape)};
        cbfp::CudaColumnBlockMatrix gpu(rows, 3, surveys);
        submit(gpu, v);
        gpu.synchronize();  // must not throw
        bool clean = true;
        for (std::size_t j = 0; j < 3; ++j) {
            clean = clean && 0 == gpu.column_contradictions(j);
        }
        check(clean, "an honest batch reports no contradiction");
    }
}

// Reading the caller's buffer in place, with an InputRead handle standing in
// for the staging copy that used to make the lifetime question go away.
void
test_in_place_input()
{
    const std::size_t rows = 2048, cols = 16, n = rows * cols;
    std::vector<double> host(n);
    for (auto &x : host) x = random_value(10);

    // The reference: the same batches through the staging path.
    cbfp::ColumnBlockMatrix cpu(rows, cols);
    submit(cpu, host);
    submit(cpu, host);

    cbfp::CudaColumnBlockMatrix staged(rows, cols);
    staged.reserve_like(cpu);
    submit(staged, host);
    submit(staged, host);

    double *pinned = nullptr;
    check(cudaSuccess == cudaHostAlloc(&pinned, n * sizeof(double),
                                       cudaHostAllocMapped),
          "pinned mapped allocation succeeds");
    std::memcpy(pinned, host.data(), n * sizeof(double));

    cbfp::CudaColumnBlockMatrix in_place(rows, cols);
    in_place.reserve_like(cpu);
    cbfp::InputRead r1 = in_place.add_matrix_col_major(pinned);
    // The buffer must not be rewritten until the handle clears; wait, then
    // resubmit the same contents so the two accumulators see the same batches.
    r1.wait();
    check(r1.ready(), "the handle reports ready once waited on");
    cbfp::InputRead r2 = in_place.add_matrix_col_major(pinned);
    r2.wait();

    std::size_t mismatches = 0;
    for (std::size_t j = 0; j < cols; ++j) {
        if (in_place.column_exponent(j) != staged.column_exponent(j) ||
            in_place.column_limbs(j) != staged.column_limbs(j)) {
            ++mismatches;
            continue;
        }
        const std::vector<std::uint64_t> a = staged.download_column(j);
        const std::vector<std::uint64_t> b = in_place.download_column(j);
        if (a != b) ++mismatches;
    }
    check(mismatches == 0, "in-place matches the staged path: " +
                               std::to_string(mismatches) + " mismatches");

    // A default-constructed handle is already clear, and so is one from an
    // empty accumulator.
    {
        cbfp::InputRead none;
        check(none.ready(), "a default handle is ready");
        none.wait();
    }

    // Move semantics: the moved-from handle must not double-destroy its event.
    {
        cbfp::InputRead a = in_place.add_matrix_col_major(pinned);
        cbfp::InputRead b = std::move(a);
        check(a.ready(), "a moved-from handle is inert");
        b.wait();
        cbfp::InputRead c;
        c = std::move(b);
        c.wait();
        check(c.ready(), "move assignment carries the event");
    }

    // Pageable memory is refused rather than left to fault in the kernel.
    {
        cbfp::CudaColumnBlockMatrix g(rows, cols);
        g.reserve_like(cpu);
        bool threw = false;
        try {
            // Deliberately NOT through submit(), which would pin it first --
            // this is the ordinary vector a caller might reach for.
            g.add_matrix_col_major(host.data());
        } catch (const std::invalid_argument &) {
            threw = true;
        }
        check(threw, "pageable input is rejected, not copied");

        // And the accumulator is still usable afterwards, so the check is a
        // rejection and not a wound.
        submit(g, host);
        g.synchronize();
        check(g.column_limbs(0) >= 1, "the accumulator survives the rejection");
    }

    // Registered rather than allocated: cudaHostRegister is the other way to
    // get memory this path accepts.
    {
        std::vector<double> own(n);
        for (auto &x : own) x = random_value(6);
        if (cudaSuccess == cudaHostRegister(own.data(), n * sizeof(double),
                                            cudaHostRegisterMapped)) {
            cbfp::ColumnBlockMatrix c2(rows, cols);
            submit(c2, own);
            cbfp::CudaColumnBlockMatrix g(rows, cols);
            g.reserve_like(c2);
            cbfp::InputRead r = g.add_matrix_col_major(own.data());
            r.wait();
            std::size_t bad = 0;
            for (std::size_t j = 0; j < cols; ++j) {
                const std::vector<std::uint64_t> dev = g.download_column(j);
                for (std::size_t i = 0; i < rows; ++i) {
                    const std::vector<std::uint64_t> h = c2.entry_limbs(i, j);
                    for (std::size_t k = 0; k < h.size(); ++k)
                        if (h[k] != dev[k * rows + i]) ++bad;
                }
            }
            check(bad == 0, "cudaHostRegister'd memory reads in place: " +
                                std::to_string(bad) + " mismatches");
            cudaHostUnregister(own.data());
        }
    }

    cudaFreeHost(pinned);
}

// Folding several device-resident matrices into one pass changes only the
// traffic. Exact accumulation does not care about grouping, so the stored
// limbs must be identical to submitting them one at a time.
void
test_device_fold()
{
    for (int variant = 0; variant < 3; ++variant) {
        const std::size_t rows = variant == 2 ? 257 : 4096;
        const std::size_t cols = variant == 2 ? 9 : 32;
        const std::size_t count = variant == 1 ? 8 : 3;
        const int spread = variant == 1 ? 300 : 8;

        std::vector<std::vector<double>> batches(
            count, std::vector<double>(rows * cols));
        for (auto &v : batches)
            for (auto &x : v) x = random_value(spread);

        // Sequential reference, one matrix per launch.
        cbfp::ColumnBlockMatrix cpu(rows, cols);
        for (auto &v : batches) submit(cpu, v);
        cbfp::CudaColumnBlockMatrix seq(rows, cols);
        seq.reserve_like(cpu);

        std::vector<double *> dev(count);
        std::vector<const double *> devc(count);
        for (std::size_t k = 0; k < count; ++k) {
            cudaMalloc(&dev[k], rows * cols * sizeof(double));
            cudaMemcpy(dev[k], batches[k].data(), rows * cols * sizeof(double),
                       cudaMemcpyHostToDevice);
            devc[k] = dev[k];
            seq.add_matrix_col_major_device(dev[k]);
        }
        seq.synchronize();

        cbfp::CudaColumnBlockMatrix fold(rows, cols);
        fold.reserve_like(cpu);
        fold.add_matrices_col_major_device(devc.data(), count);
        fold.synchronize();

        std::size_t bad = 0;
        for (std::size_t j = 0; j < cols; ++j) {
            if (fold.column_exponent(j) != seq.column_exponent(j) ||
                fold.column_limbs(j) != seq.column_limbs(j)) {
                ++bad;
                continue;
            }
            if (fold.download_column(j) != seq.download_column(j)) ++bad;
        }
        check(bad == 0, "device fold matches sequential, variant " +
                            std::to_string(variant) + ": " +
                            std::to_string(bad) + " mismatches");
        for (auto *p : dev) cudaFree(p);
    }

    // An adaptively sized accumulator has to survey every matrix and fit once
    // before adding any of them, or a widen partway through would leave the
    // earlier ones written into a column of the wrong shape.
    {
        const std::size_t rows = 512, cols = 4, count = 3;
        std::vector<std::vector<double>> batches(
            count, std::vector<double>(rows * cols));
        for (std::size_t k = 0; k < count; ++k)
            for (auto &x : batches[k])
                x = std::ldexp(1.0 + (double)(rng() % 100), -40 * (int)k);

        cbfp::ColumnBlockMatrix cpu(rows, cols);
        for (auto &v : batches) submit(cpu, v);

        std::vector<double *> dev(count);
        std::vector<const double *> devc(count);
        for (std::size_t k = 0; k < count; ++k) {
            cudaMalloc(&dev[k], rows * cols * sizeof(double));
            cudaMemcpy(dev[k], batches[k].data(), rows * cols * sizeof(double),
                       cudaMemcpyHostToDevice);
            devc[k] = dev[k];
        }
        // Reserved wide enough only for the first matrix, so the fold must
        // widen and lower the exponent from the aggregate of all three.
        cbfp::CudaColumnBlockMatrix g(rows, cols);
        for (std::size_t j = 0; j < cols; ++j) g.reserve_column(j, 0, 64);
        g.add_matrices_col_major_device(devc.data(), count);
        g.synchronize();
        check(g.column_exponent(0) < -60,
              "an adaptive fold lowers the exponent to the aggregate");
        std::size_t bad = 0;
        for (std::size_t j = 0; j < cols; ++j) {
            const std::vector<std::uint64_t> d = g.download_column(j);
            for (std::size_t i = 0; i < rows; ++i) {
                const std::vector<std::uint64_t> h = cpu.entry_limbs(i, j);
                if (g.column_exponent(j) != cpu.column_exponent(j) ||
                    g.column_limbs(j) != h.size()) { ++bad; break; }
                for (std::size_t k = 0; k < h.size(); ++k)
                    if (h[k] != d[k * rows + i]) ++bad;
            }
        }
        check(bad == 0, "and gets the same answer: " + std::to_string(bad) +
                            " mismatches");
        for (auto *p : dev) cudaFree(p);
    }

    // More than the parameter block holds is an error, not a silent truncation.
    {
        cbfp::ColumnBlockMatrix cpu(64, 2);
        std::vector<double> v(64 * 2, 1.0);
        submit(cpu, v);
        cbfp::CudaColumnBlockMatrix g(64, 2);
        g.reserve_like(cpu);
        double *d = nullptr;
        cudaMalloc(&d, 64 * 2 * sizeof(double));
        std::vector<const double *> many(17, d);
        bool threw = false;
        try {
            g.add_matrices_col_major_device(many.data(), many.size());
        } catch (const std::invalid_argument &) {
            threw = true;
        }
        check(threw, "folding more matrices than fit is rejected");
        cudaFree(d);
    }
}

// The pattern the device path exists for: a producer cycling device buffers,
// filling one while the accumulator reads another. The handle is what makes it
// possible to wait for the one buffer about to be overwritten instead of
// draining the whole stream.
void
test_device_ping_pong()
{
    const std::size_t rows = 4096, cols = 16, n = rows * cols;
    const int batches = 6;

    std::vector<std::vector<double>> host(batches, std::vector<double>(n));
    for (auto &v : host)
        for (auto &x : v) x = random_value(10);

    cbfp::ColumnBlockMatrix cpu(rows, cols);
    for (auto &v : host) submit(cpu, v);

    // Two device buffers in rotation, which is all a real producer would keep.
    double *buf[2] = {nullptr, nullptr};
    for (int i = 0; i < 2; ++i) cudaMalloc(&buf[i], n * sizeof(double));

    cbfp::CudaColumnBlockMatrix g(rows, cols);
    g.reserve_like(cpu);
    cbfp::InputRead h[2];
    for (int r = 0; r < batches; ++r) {
        const int slot = r & 1;
        // Wait only for the buffer about to be refilled -- never for the
        // stream. On the first two rounds the handles are empty and this is a
        // no-op, which is what a default-constructed handle is for.
        h[slot].wait();
        cudaMemcpy(buf[slot], host[r].data(), n * sizeof(double),
                   cudaMemcpyHostToDevice);
        h[slot] = g.add_matrix_col_major_device(buf[slot]);
    }
    h[0].wait();
    h[1].wait();
    g.synchronize();

    std::size_t bad = 0;
    for (std::size_t j = 0; j < cols; ++j) {
        if (g.column_exponent(j) != cpu.column_exponent(j) ||
            g.column_limbs(j) != cpu.column_limbs(j)) {
            ++bad;
            continue;
        }
        const std::vector<std::uint64_t> dev = g.download_column(j);
        for (std::size_t i = 0; i < rows; ++i) {
            const std::vector<std::uint64_t> hl = cpu.entry_limbs(i, j);
            for (std::size_t k = 0; k < hl.size(); ++k)
                if (hl[k] != dev[k * rows + i]) ++bad;
        }
    }
    check(bad == 0, "two device buffers in rotation match sequential: " +
                        std::to_string(bad) + " mismatches");

    // The handle really does gate the read: a buffer whose handle has cleared
    // may be overwritten, and doing so must not disturb what was accumulated.
    {
        const std::vector<std::uint64_t> before = g.download_column(0);
        h[0].wait();
        h[1].wait();
        cudaMemset(buf[0], 0xFF, n * sizeof(double));
        cudaMemset(buf[1], 0xFF, n * sizeof(double));
        cudaDeviceSynchronize();
        check(g.download_column(0) == before,
              "overwriting a cleared buffer does not disturb the accumulator");
    }

    for (int i = 0; i < 2; ++i) cudaFree(buf[i]);
}

// Zeroing a freshly reserved column is deferred to the first launch or
// readback, so that a whole pre-sized accumulator is zeroed in one kernel
// rather than one cudaMemsetAsync per limb position. A reader must not be able
// to see that: the limbs read as zero whether or not anything has been added.
void
test_deferred_zeroing()
{
    const std::size_t rows = 300, cols = 5;

    // Read back before anything is accumulated -- the readback path has to
    // settle the pending zeroing itself.
    {
        cbfp::CudaColumnBlockMatrix g(rows, cols);
        for (std::size_t j = 0; j < cols; ++j) g.reserve_column(j, 0, 256);
        bool zero = true;
        for (std::size_t j = 0; j < cols; ++j) {
            const std::vector<std::uint64_t> col = g.download_column(j);
            for (std::uint64_t w : col) zero = zero && (0 == w);
            for (std::size_t i = 0; i < rows; i += 97) {
                for (std::uint64_t w : g.entry_limbs(i, j)) {
                    zero = zero && (0 == w);
                }
            }
        }
        check(zero, "a reserved column reads as zero before any accumulation");
    }

    // And a column reserved but never written stays zero while its neighbours
    // are accumulated into, which is what would break if the flush zeroed only
    // what the launch happens to touch.
    {
        std::vector<double> v(rows * cols, 0.0);
        for (std::size_t j = 0; j < cols; ++j) {
            for (std::size_t i = 0; i < rows; ++i) {
                v[j * rows + i] = (j == 2) ? 0.0 : random_value(6);
            }
        }
        cbfp::ColumnBlockMatrix cpu(rows, cols);
        submit(cpu, v);
        cbfp::CudaColumnBlockMatrix g(rows, cols);
        g.reserve_like(cpu);
        submit(g, v);
        g.synchronize();

        bool zero = true;
        for (std::uint64_t w : g.download_column(2)) zero = zero && (0 == w);
        check(zero, "an all-zero column stays zero through a launch");

        std::size_t bad = 0;
        for (std::size_t j = 0; j < cols; ++j) {
            const std::vector<std::uint64_t> dev = g.download_column(j);
            for (std::size_t i = 0; i < rows; ++i) {
                const std::vector<std::uint64_t> h = cpu.entry_limbs(i, j);
                for (std::size_t k = 0; k < h.size(); ++k) {
                    if (h[k] != dev[k * rows + i]) ++bad;
                }
            }
        }
        check(bad == 0, "and the rest matches the host: " +
                            std::to_string(bad) + " mismatches");
    }
}

// The device survey has to agree with the host one field for field, or a
// producer computing it on whichever side its data happens to be on would
// size an accumulator differently depending on where the survey ran.
void
test_device_survey()
{
    const std::size_t rows = 3000, cols = 17, n = rows * cols;
    std::vector<double> src(n);
    for (std::size_t j = 0; j < cols; ++j) {
        for (std::size_t i = 0; i < rows; ++i) {
            double v = random_value(120);
            if (1 == j) v = 0.0;                       // an all-zero column
            if (2 == j && 0 == (i % 7)) v = 0.0;       // zeros mixed in
            if (3 == j) v = std::ldexp(1.0, -1074);    // all subnormal
            src[j * rows + i] = v;
        }
    }
    double *dev = nullptr;
    cudaMalloc(&dev, n * sizeof(double));
    cudaMemcpy(dev, src.data(), n * sizeof(double), cudaMemcpyHostToDevice);

    std::vector<cbfp::Survey> host(cols), device(cols);
    cbfp::survey_matrix_col_major(host.data(), src.data(), rows, cols);

    // The device writes `out`, so it goes to device memory and comes back by
    // an explicit copy the caller chose to make.
    cbfp::Survey *d_out = nullptr;
    cudaMalloc(&d_out, cols * sizeof(cbfp::Survey));
    cbfp::survey_matrix_col_major_device(d_out, dev, rows, cols);
    cudaMemcpy(device.data(), d_out, cols * sizeof(cbfp::Survey),
               cudaMemcpyDeviceToHost);
    cudaFree(d_out);

    std::size_t bad = 0;
    for (std::size_t j = 0; j < cols; ++j) {
        if (host[j].min_exponent != device[j].min_exponent ||
            host[j].max_top != device[j].max_top ||
            host[j].any != device[j].any ||
            host[j].nonfinite != device[j].nonfinite) {
            ++bad;
        }
    }
    check(bad == 0, "device survey matches the host one: " +
                        std::to_string(bad) + " of " + std::to_string(cols) +
                        " columns differ");
    check(!host[1].any && !device[1].any, "and both report the empty column");

    // Ordinary pageable output is refused rather than faulting in the kernel.
    {
        std::vector<cbfp::Survey> pageable(cols);
        bool threw = false;
        try {
            cbfp::survey_matrix_col_major_device(pageable.data(), dev, rows,
                                                 cols);
        } catch (const std::invalid_argument &) {
            threw = true;
        }
        check(threw, "a pageable output pointer is rejected");
    }

    // A survey taken on the device must size an accumulator the same way one
    // taken on the host does -- which is the whole point of sending it ahead.
    {
        std::vector<std::vector<cbfp::Survey>> from_dev{device};
        std::vector<std::vector<cbfp::Survey>> from_host{host};
        cbfp::ColumnBlockMatrix a(rows, cols, from_dev);
        cbfp::ColumnBlockMatrix b(rows, cols, from_host);
        bool same = true;
        for (std::size_t j = 0; j < cols; ++j) {
            same = same && a.column_exponent(j) == b.column_exponent(j) &&
                   a.column_limbs(j) == b.column_limbs(j);
        }
        check(same, "and reserves identically to it");
    }

    // Non-finite input is reported, not silently surveyed.
    {
        std::vector<double> bad_in(rows * 2, 1.0);
        bad_in[rows + 5] = std::numeric_limits<double>::infinity();
        double *d2 = nullptr;
        cudaMalloc(&d2, rows * 2 * sizeof(double));
        cudaMemcpy(d2, bad_in.data(), rows * 2 * sizeof(double),
                   cudaMemcpyHostToDevice);
        // Mapped pinned output: the other kind of memory the device can
        // write, and how a caller asks for the answer on the host.
        cbfp::Survey *mapped = nullptr;
        cudaHostAlloc(&mapped, 2 * sizeof(cbfp::Survey), cudaHostAllocMapped);
        cbfp::survey_matrix_col_major_device(mapped, d2, rows, 2);
        check(!mapped[0].nonfinite && mapped[1].nonfinite,
              "the device survey flags a non-finite column and only that one");
        check(mapped[0].any && 0 == mapped[1].min_exponent,
              "and a mapped pinned output is readable on return");
        cudaFreeHost(mapped);
        cudaFree(d2);
    }
    cudaFree(dev);
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
    test_rescale();
    test_acquired_input();
    test_errors();
    test_symmetric();
    test_presized();
    test_device_detector();
    test_in_place_input();
    test_device_fold();
    test_device_ping_pong();
    test_deferred_zeroing();
    test_device_survey();

    std::printf("%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
