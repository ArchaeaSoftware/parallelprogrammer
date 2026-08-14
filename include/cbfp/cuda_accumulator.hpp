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
// Both column parameters adapt, as they do on the CPU, but between launches
// rather than inside one. A column that needs more limbs appends them and
// sign-fills them, with nothing already allocated moved; a column that needs a
// lower exponent has every entry shifted left to match. Neither can happen
// from inside a kernel, so both are driven from the host once the survey has
// said what the pending batch needs.
//
// reserve_for is therefore an optimization here exactly as it is on the CPU:
// results are the same without it, it just avoids the growing and shifting.
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "cbfp/limbs.hpp"
#include "cbfp/survey.hpp"
#include "cbfp/uplo.hpp"

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

// Signals when the device has finished reading an input buffer that was
// submitted to be read in place, which is what tells a caller its buffer is
// writable again.
//
// The hazard this exists for is real and was measured: the kernel streams the
// buffer over PCIe and is still reading it after the submitting call returns,
// so a caller refilling it corrupts about 63% of entries. An earlier design
// tried to infer safety from the pointer's memory type, which gave one
// function two lifetime contracts and nothing in the signature to tell them
// apart. Handing back a handle puts the obligation where the caller can see it.
//
// Move-only, and owns a CUDA event. Destroying one without waiting is not an
// error -- it is a statement that the buffer will not be written again.
class InputRead {
public:
    InputRead() = default;
    ~InputRead();
    InputRead(InputRead &&other) noexcept;
    InputRead &operator=(InputRead &&other) noexcept;
    InputRead(const InputRead &) = delete;
    InputRead &operator=(const InputRead &) = delete;

    // Blocks until the device has finished reading. A default-constructed
    // handle, or one whose submission was empty, returns at once.
    void wait();

    // The same question without blocking.
    bool ready() const;

private:
    friend class CudaColumnBlockMatrix;
    explicit InputRead(CUevent_st *ev) : ev_(ev) {}
    CUevent_st *ev_ = nullptr;
};

class CudaColumnBlockMatrix {
public:
    CudaColumnBlockMatrix(std::size_t rows, std::size_t cols);

    // Symmetric n x n, storing one triangle, exactly as ColumnBlockMatrix
    // does. Column j holds n-j entries (Lower) or j+1 (Upper), and each
    // column's stored rows stay contiguous, so a warp still reads 32
    // consecutive rows of one limb column as a single coalesced access.
    //
    // The kernels read only the stored slice of an input matrix, so the
    // bytes crossing PCIe halve along with the device memory. The input is
    // taken to be symmetric and that is not checked.
    CudaColumnBlockMatrix(std::size_t n, Uplo uplo);

    // Pre-sized from the surveys of every matrix that will be accumulated,
    // exactly as ColumnBlockMatrix is: `surveys` is indexed by matrix and
    // then by column, so its outer size is how many will arrive.
    //
    // This buys more here than on the CPU. The device path's survey is a
    // kernel whose verdict has to reach the host before the accumulate may
    // launch -- a drain that costs 19.9 us a batch however small the batch
    // is. Told the extents in advance there is nothing to ask, so the survey
    // kernel, its round trip and that fixed cost all go away, and
    // add_matrix_col_major_device stops synchronizing altogether.
    //
    // The metadata is taken on trust. A batch that contradicts it is
    // detected by the accumulate kernel rather than prevented, and reported
    // at the next synchronization -- see column_contradictions().
    CudaColumnBlockMatrix(std::size_t rows, std::size_t cols,
                          const std::vector<std::vector<Survey>> &surveys);

    CudaColumnBlockMatrix(std::size_t n, Uplo uplo,
                          const std::vector<std::vector<Survey>> &surveys);

    // Declared, not implicit: the worker pool is held by unique_ptr to an
    // incomplete type, so the destructor must be defined where that is.
    ~CudaColumnBlockMatrix();

    CudaColumnBlockMatrix(const CudaColumnBlockMatrix &) = delete;
    CudaColumnBlockMatrix &operator=(const CudaColumnBlockMatrix &) = delete;

    std::size_t rows() const { return rows_; }

    std::size_t cols() const { return cols_; }

    bool symmetric() const { return symmetric_; }

    // Only meaningful when symmetric().
    Uplo uplo() const { return uplo_; }

    // Entries stored in column j, and the logical row the first of them is.
    // rows() and 0 unless symmetric.
    std::size_t column_rows(std::size_t j) const;
    std::size_t column_first_row(std::size_t j) const;

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

    // --- input memory ------------------------------------------------------

    // Borrows the next input buffer, sized rows() * cols() doubles. Fill it
    // column-major and pass it to add_matrix_col_major, which recognises it
    // and lets the kernel stream it over PCIe in place -- no staging copy, and
    // no device-side copy of the input at all.
    //
    // Blocks only if the device is still reading the buffer being handed back,
    // so with two in rotation the host fills one while the device streams the
    // other. That rotation is the point: reading in place means the kernel is
    // still using a buffer after the call that submitted it returned, and a
    // caller refilling it would be writing under the device. Measured, doing
    // so corrupts about 63% of entries.
    //
    // The buffer belongs to the accumulator and stays valid until the next
    // acquire_input. Ordinary host memory still works everywhere and is staged
    // through the same buffers; this only removes that copy.
    double *acquire_input();

    // --- accumulation ------------------------------------------------------

    // A += B, where B is column-major on the host: column j begins at
    // b + j*col_stride and its rows are contiguous (0 means tightly packed).
    // Throws std::domain_error on inf/NaN, and std::runtime_error if a column
    // was reserved too narrow or at too high an exponent for these values.
    void add_matrix_col_major(const double *b, std::size_t col_stride = 0);

    // A += B, reading B where it lies instead of staging a copy of it.
    //
    // B must be page-locked and device-mapped -- cudaHostAlloc with
    // cudaHostAllocMapped, or cudaHostRegister with cudaHostRegisterMapped --
    // because the kernel dereferences it on the device. That is checked here
    // rather than left to fault at the launch.
    //
    // B is still being read when this returns. Do not write to it until the
    // returned handle says the read has finished. This is the same path
    // acquire_input offers, for a caller that would rather own the buffer; the
    // staging copy it avoids costs 1331 us at 65536x64, which is most of what
    // separates the two rows of the table above.
    InputRead add_matrix_col_major_in_place(const double *b,
                                            std::size_t col_stride = 0);

    // The same, for input already resident in device memory.
    void add_matrix_col_major_device(const double *b,
                                     std::size_t col_stride = 0);

    // Blocks until every submitted accumulation has finished. Accumulation is
    // asynchronous: a call returns once the work is queued, so a caller
    // streaming batches keeps the device busy without doing anything special.
    // Readback synchronizes implicitly, so this is only needed for timing or
    // before reusing a caller-owned device buffer.
    void synchronize() const;

    // What the accumulate kernel found that contradicted a column's sizing,
    // as kernels.hpp's kBad* bits, or 0. Accumulated across every batch, the
    // way occupancy is. Reading it synchronizes.
    //
    // synchronize() throws on a nonzero value, so a caller normally learns
    // about a contradiction rather than asking; this is here for a caller
    // that wants to know which column without catching.
    unsigned column_contradictions(std::size_t j) const;

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

        // Entries stored in this column and the logical row of the first,
        // held per column because a triangular column has its own length and
        // every allocation, widen and rescale is sized from it.
        std::size_t rows = 0;
        std::size_t first_row = 0;
        // One device allocation per limb position, mirroring the CPU's
        // vector<LimbColumn>. Widening is then an append: the existing
        // allocations are not touched at all.
        std::vector<limbs::limb_t *> bases;
        limbs::limb_t **dev_bases = nullptr;  // the same array, device-side

        // The same derived width bound ColumnBlockMatrix keeps, and for the
        // same reason: no entry can exceed count * 2^max_addend_bits, so that
        // bound plus a sign bit says how wide the column must be. Checking
        // only the incoming addend is not enough -- an addend does not grow
        // with the batch count but the sum does, and a column accumulated past
        // what it was reserved for wraps silently.
        std::size_t max_addend_bits = 0;
        std::size_t add_count = 0;

        // Highest limb position the accumulate has ever disturbed, reported
        // by the kernel rather than derived. -1 until something is added.
        // Monotonic: cancellation can shrink the value but not this, so it
        // stays an upper bound on how much of the column is in use.
        int max_limb_used = -1;
    };

    // Two staging slots, so the host can prepare batch N+1 while the device is
    // still accumulating batch N.
    //
    // The buffer is mapped rather than merely pinned, and there is no device
    // copy of it: the kernel reads host memory directly. Measured at 65536x64,
    // copying 33.6 MB and then reading it from device memory takes 1768 us,
    // while reading it in place takes 1257 -- exactly what the transfer alone
    // costs, so the accumulator work hides entirely behind the bus. That only
    // holds because the input is read once: the host survey is what means the
    // device never needs a second look at it.
    struct Slot {
        double *mapped = nullptr;
        CUevent_st *ev_done = nullptr;  // last kernel to read it finished
        bool in_flight = false;
    };

    void check_index(std::size_t i, std::size_t j) const;
    void reserve_from_surveys(const std::vector<std::vector<Survey>> &s);
    // Throws if any column reported a contradiction. Const because it only
    // reads mapped memory the device wrote.
    void report_contradictions() const;
    // Logical (i, j) to the column that stores it and the slot within it.
    void locate(std::size_t i, std::size_t j, std::size_t &col,
                std::size_t &slot) const;
    void init_columns();
    void accumulate_device(const double *b, std::size_t col_stride);
    void launch_accumulate(const double *b, std::size_t col_stride);
    void validate_host_survey(const double *b, std::size_t col_stride);
    void sync_descriptors();
    void harvest_occupancy();
    void reserve_from_extents(const int *low, const int *high,
                              const char *any, std::size_t count);
    void require_fit(std::size_t j, int min_exponent, int max_top);
    void grow_column(std::size_t j, std::size_t needed);
    void rescale_column(std::size_t j, int new_exponent);
    void ensure_slots(std::size_t words);

    std::size_t rows_;
    std::size_t cols_;
    bool symmetric_ = false;
    Uplo uplo_ = Uplo::Lower;

    // Deliberately no skew between limb columns, unlike the CPU. The CPU needs
    // it because eight lanes touch nlimbs arrays at the same row offset in
    // succession and collide in an 8-way L1 set. A warp instead reads 32
    // consecutive rows of one limb column as a single coalesced transaction,
    // and visits limb positions sequentially within a thread, so there is no
    // equivalent collision -- and the allocator's 256-byte alignment, which is
    // what coalescing actually needs, comes for free.
    std::vector<Column> cols_state_;

    // One stream. There is no longer a transfer to overlap with compute --
    // the kernel does the transfer, by reading host memory as it goes.
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

    // The same pair for the detector's verdict, carried across by the same
    // elected block in the same pass.
    unsigned *flags_device_ = nullptr;
    unsigned *flags_host_ = nullptr;

    void *descriptors_ = nullptr;  // device ColumnDesc[]
    bool descriptors_stale_ = true;
    void *survey_out_ =
        nullptr;  // device survey results, for device-side input

    // Per-column stored length and first row, device-side. Fixed at
    // construction and never rewritten, so both kernels can read it without
    // the staleness the descriptors have to manage.
    void *shapes_ = nullptr;  // device ColumnShape[]

    // Set when the constructor was given surveys. The accumulate path then
    // skips the survey entirely -- on this path that means not launching a
    // kernel and not draining the stream to read its answer.
    bool presized_ = false;
    std::size_t declared_matrices_ = 0;
    std::size_t submitted_matrices_ = 0;

    Slot slots_[2];
    int slot_ = 0;
    std::size_t slot_words_ = 0;
};

}  // namespace cbfp
