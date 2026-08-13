// Device-resident sibling of ColumnBlockMatrix.
//
// This is deliberately not a third kernel variant behind the CPU's dispatch.
// Swapping a CPU kernel leaves the memory and the object identical, so that is
// a true backend swap; CUDA's limb arrays live in device memory, and putting
// them behind the same function pointer would disguise host/device transfers
// as ordinary calls. What the two implementations share is the *algorithm* --
// when to rescale, when to widen, how the exponent moves -- not the memory.
// See docs/simd-design.md.
//
// Pre-sizing is mandatory here, where it is only an optimization on the CPU.
// Growing a column mid-kernel would mean reallocating device memory from
// inside a launch, so a column's exponent and width are fixed before any
// accumulation touches it, and values that would not fit are an error rather
// than a silent rescale.
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "cbfp/limbs.hpp"

// CUDA's opaque handle types, forward-declared rather than pulled in from
// cuda_runtime.h so this header stays usable from a plain C++ translation unit
// with no CUDA toolkit headers on its include path -- only the link needs
// cudart. cudaStream_t and cudaEvent_t are exactly these pointer types, and
// cuda_accumulator.cu static_asserts as much, so drift would be a compile
// error rather than a silent one.
struct CUstream_st;
struct CUevent_st;

namespace cbfp {

class ColumnBlockMatrix;

// True if a CUDA device is present and usable. Everything below throws
// std::runtime_error if it is not.
bool
cuda_available();

class CudaColumnBlockMatrix {
public:
    CudaColumnBlockMatrix(std::size_t rows, std::size_t cols);
    ~CudaColumnBlockMatrix();

    CudaColumnBlockMatrix(const CudaColumnBlockMatrix &) = delete;
    CudaColumnBlockMatrix &operator=(const CudaColumnBlockMatrix &) = delete;

    std::size_t rows() const { return rows_; }

    std::size_t cols() const { return cols_; }

    // --- pre-sizing --------------------------------------------------------

    // Fix column j's scale and width, and allocate its limb arrays. `exponent`
    // is the lowest bit weight that will be needed and `bits` the total width.
    // May be called once per column, before anything is accumulated into it.
    void reserve_column(std::size_t j, int exponent, std::size_t bits);

    // Take every column's scale and width from a CPU accumulator that has
    // already seen the data. This is the natural way to drive the device path
    // and what makes the two directly comparable: given the same exponent and
    // width, both must hold bit-identical limbs.
    void reserve_like(const ColumnBlockMatrix &cpu);

    // Pre-size every column from a representative batch about to be
    // accumulated `count` times, the same arithmetic as
    // ColumnBlockMatrix::reserve_for. This is what lets the device container
    // stand on its own: reserve_like needs a CPU accumulator that has already
    // seen the data, which means doing the whole computation twice.
    //
    // `b` is column-major on the host. A column that is entirely zero still
    // gets a minimal reservation, because an unreserved column is an error
    // here rather than something that can grow on first use.
    void reserve_for(const double *b, std::size_t count = 1,
                     std::size_t col_stride = 0);

    // The same, for a representative batch already resident on the device.
    void reserve_for_device(const double *b, std::size_t count = 1,
                            std::size_t col_stride = 0);

    // --- accumulation ------------------------------------------------------

    // A += B, where B is column-major on the host: column j begins at
    // b + j*col_stride and its rows are contiguous (0 means tightly packed).
    // Throws std::domain_error on inf/NaN, and std::runtime_error if a column
    // was reserved too narrow or at too high an exponent for these values.
    void add_matrix_col_major(const double *b, std::size_t col_stride = 0);

    // The same, for input already resident in device memory.
    void add_matrix_col_major_device(const double *b,
                                     std::size_t col_stride = 0);

    // Blocks until every submitted accumulation has finished. Accumulation is
    // asynchronous: a call returns once the work is queued, so a caller
    // streaming batches keeps the device busy without doing anything special.
    // Readback synchronizes implicitly, so this is only needed for timing or
    // before reusing a caller-owned device buffer.
    void synchronize() const;

    // --- readback ----------------------------------------------------------

    int column_exponent(std::size_t j) const;

    std::size_t column_limbs(std::size_t j) const;

    // Entry (i, j)'s stored two's complement limbs, copied back from the
    // device. The counterpart of ColumnBlockMatrix::entry_limbs.
    std::vector<limbs::limb_t> entry_limbs(std::size_t i, std::size_t j) const;

    // Column j's limb arrays in bulk: out[k*rows() + i] is entry i's k-th
    // limb. One copy per limb position rather than one per entry, which is
    // what makes a whole-matrix comparison practical.
    std::vector<limbs::limb_t> download_column(std::size_t j) const;

    // Highest limb position column j has ever needed, as measured by the
    // accumulate rather than derived from a worst case. -1 if nothing has
    // been added yet. Readback synchronizes, so this reflects every batch
    // submitted so far.
    int column_occupancy(std::size_t j) const;

    std::size_t memory_bytes() const;

private:
    struct Column {
        int exponent = 0;
        bool reserved = false;
        std::size_t nlimbs = 0;
        // One device allocation per limb position, mirroring the CPU's
        // vector<LimbColumn>. Widening is then an append: the existing
        // allocations are not touched at all.
        std::vector<limbs::limb_t *> bases;
        limbs::limb_t **dev_bases = nullptr;  // the same array, device-side

        // Highest limb position the accumulate has ever disturbed, reported
        // by the kernel rather than derived. -1 until something is added.
        // Monotonic: cancellation can shrink the value but not this, so it
        // stays an upper bound on how much of the column is in use.
        int max_limb_used = -1;
    };

    // Two staging slots, so the host can prepare batch N+1 while the device is
    // still accumulating batch N. Each owns a pinned host buffer -- measured,
    // an async copy from pageable memory blocks the host for the whole
    // transfer and there is no window to survey in -- the device buffer the
    // kernel reads, and an event marking when the last kernel to read that
    // buffer finished, so a slot is only reused once it is genuinely free.
    struct Slot {
        double *pinned = nullptr;
        double *device = nullptr;
        CUevent_st *ev_copied = nullptr;  // H2D into `device` finished
        CUevent_st *ev_done = nullptr;    // last kernel to read it finished
        bool in_flight = false;
    };

    void check_index(std::size_t i, std::size_t j) const;
    void accumulate_device(const double *b, std::size_t col_stride);
    void launch_accumulate(const double *b, std::size_t col_stride);
    void validate_host_survey(const double *packed);
    void sync_descriptors();
    void harvest_occupancy();
    void reserve_from_extents(const long long *low, const long long *high,
                              const char *any, std::size_t count);
    void ensure_slots(std::size_t words);

    std::size_t rows_;
    std::size_t cols_;

    // Deliberately no skew between limb columns, unlike the CPU. The CPU needs
    // it because eight lanes touch nlimbs arrays at the same row offset in
    // succession and collide in an 8-way L1 set. A warp instead reads 32
    // consecutive rows of one limb column as a single coalesced transaction,
    // and visits limb positions sequentially within a thread, so there is no
    // equivalent collision -- and the allocator's 256-byte alignment, which is
    // what coalescing actually needs, comes for free.
    std::vector<Column> cols_state_;

    // Copy and compute are separate streams so batch N+1's transfer runs on
    // the copy engine while batch N is still accumulating. On one stream they
    // serialize, and the transfer is the longer of the two.
    CUstream_st *st_copy_ = nullptr;
    CUstream_st *st_compute_ = nullptr;
    // These two stay void*: they point at types defined inside the .cu, which
    // is where they belong -- the descriptor layout is not this header's
    // business.
    // The accumulate's per-column results come back through mapped host
    // memory: every block reduces into occupancy_device_ with ordinary device
    // atomics, then the last block to finish copies the finished array across
    // and re-arms it. Mapped rather than copied because it is one value per
    // column -- measured, that wins below ~64 values and loses badly above.
    int *occupancy_device_ = nullptr;  // staging, device
    int *occupancy_host_ = nullptr;    // mapped, written by the last block
    unsigned *ticket_ = nullptr;       // device, elects that block

    void *descriptors_ = nullptr;  // device ColumnDesc[]
    bool descriptors_stale_ = true;
    void *survey_out_ =
        nullptr;  // device survey results, for device-side input

    Slot slots_[2];
    int slot_ = 0;
    std::size_t slot_words_ = 0;
};

}  // namespace cbfp
