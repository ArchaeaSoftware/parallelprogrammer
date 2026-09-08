// Device-resident sibling of AccumulationMatrix.
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

#include "truesum/limbs.hpp"
#include "truesum/survey.hpp"
#include "truesum/uplo.hpp"

// CUDA's opaque handle types, forward-declared rather than pulled in from
// cuda_runtime.h so this header stays usable from a plain C++ translation unit
// with no CUDA toolkit headers on its include path -- only the link needs
// cudart. cudaStream_t and cudaEvent_t are exactly these pointer types, and
// cuda_accumulation_matrix.cu static_asserts as much, so drift would be a
// compile error rather than a silent one.
struct CUstream_st;
struct CUevent_st;

namespace truesum {

class AccumulationMatrix;

// True if a CUDA device is present and usable. Everything below throws
// std::runtime_error if it is not.
bool
cuda_available();

// Surveys a column-major matrix the device can read in place, writing `cols`
// entries. Column j begins at b + j*col_stride and its rows are
// contiguous; 0 means tightly packed.
//
// The device counterpart of survey_matrix_col_major, and the piece a producer
// needs when its data never passes through host memory -- delivered by
// GPUDirect, or computed on the device. The result is 12 bytes a column
// against 8 a matrix element, so it can be sent ahead to whoever will do the
// accumulating while the matrix itself is still in flight: 768 bytes against
// 33.6 MB at 65536x64.
//
// **Both pointers are dereferenced by the device**, so each must be device
// memory or host memory that is page-locked and mapped -- cudaMalloc, or
// cudaHostAlloc / cudaHostRegister with the mapped flag. Both are checked, and
// pageable memory is rejected rather than staged. That is the same rule the
// accumulate path follows, and it is what keeps this function free of
// allocations: it copies nothing, allocates nothing, and leaves the caller's
// allocator undisturbed between their own launches.
//
// A caller wanting the answer on the host passes a mapped pinned `out` and says
// so, rather than this function copying it there and assuming that is what was
// wanted; a caller feeding it to something else on the device pays for no copy
// at all.
//
// Asynchronous with respect to the host, as every kernel-based call here is.
// The launches go on `stream`, or on the default stream when that is null, and
// this returns as soon as they are queued. `out` holds the survey once the
// caller has synchronized that stream or waited on an event they recorded on
// it; nothing is allocated and nothing is copied, so no hidden temporary's
// lifetime forces a wait that the caller did not ask for.
void
survey_matrix_col_major_device(Survey *out, const double *b,
                               std::size_t rows, std::size_t cols,
                               std::size_t col_stride = 0,
                               CUstream_st *stream = nullptr);

class CudaAccumulationMatrix {
public:
    // `stream` is optional. Left null, the accumulation matrix creates a
    // stream of its own and destroys it in the destructor; that one is a
    // blocking stream, so a caller who queues work on the legacy null stream is
    // implicitly ordered against the accumulates and cannot be caught out by
    // it.
    //
    // Supplied, every launch goes on the caller's stream instead, and the
    // caller keeps ownership: it is not destroyed with the accumulation matrix,
    // and its flags are the caller's choice. This is how to get the overlap the
    // default gives up. A producer that queues its fills on the same stream is
    // ordered against the accumulates in both directions, needs no events, and
    // never touches the null stream -- see stream().
    CudaAccumulationMatrix(std::size_t rows, std::size_t cols,
                           CUstream_st *stream = nullptr);

    // Symmetric n x n, storing one triangle, exactly as AccumulationMatrix
    // does. Column j holds n-j entries (Lower) or j+1 (Upper), and each
    // column's stored rows stay contiguous, so a warp still reads 32
    // consecutive rows of one limb column as a single coalesced access.
    //
    // The kernels read only the stored slice of an input matrix, so the
    // bytes crossing PCIe halve along with the device memory. The input is
    // taken to be symmetric and that is not checked.
    CudaAccumulationMatrix(std::size_t n, Uplo uplo,
                           CUstream_st *stream = nullptr);

    // Pre-sized from the surveys of every matrix that will be accumulated,
    // exactly as AccumulationMatrix is: `surveys` is indexed by matrix and
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
    CudaAccumulationMatrix(std::size_t rows, std::size_t cols,
                           const std::vector<std::vector<Survey>> &surveys,
                           CUstream_st *stream = nullptr);

    CudaAccumulationMatrix(std::size_t n, Uplo uplo,
                           const std::vector<std::vector<Survey>> &surveys,
                           CUstream_st *stream = nullptr);

    // Declared, not implicit: the worker pool is held by unique_ptr to an
    // incomplete type, so the destructor must be defined where that is.
    ~CudaAccumulationMatrix();

    CudaAccumulationMatrix(const CudaAccumulationMatrix &) = delete;
    CudaAccumulationMatrix &operator=(const CudaAccumulationMatrix &) = delete;

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

    // Set column j's scale and width, and allocate its limb arrays. `exponent`
    // is the lowest bit weight that will be needed and `bits` the total width.
    // May be called once per column, before anything is accumulated into it.
    void reserve_column(std::size_t j, int exponent, std::size_t bits);

    // Take every column's scale and width from a CPU accumulation matrix that
    // has already seen the data. This is the natural way to drive the device
    // path and what makes the two directly comparable: given the same exponent
    // and width, both must hold bit-identical limbs.
    void reserve_like(const AccumulationMatrix &cpu);

    // Reserve room for matrices described by `surveys`, indexed by matrix and
    // then by column, without pre-sizing. The accumulation matrix goes on
    // surveying every matrix it is handed; what this removes is the rescale a
    // later low value would otherwise force, which on this target is an error
    // rather than an adjustment once a launch is in flight.
    void reserve_for_surveys(const std::vector<std::vector<Survey>> &surveys);

    // Pre-size every column from a representative batch about to be
    // accumulated `count` times, the same arithmetic as
    // AccumulationMatrix::reserve_for. This is what lets the device container
    // stand on its own: reserve_like needs a CPU accumulation matrix that has
    // already seen the data, which means doing the whole computation twice.
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
    // column-major and pass it to add_matrix_col_major, which recognizes it
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
    // The buffer belongs to the accumulation matrix and stays valid until the
    // next acquire_input. Ordinary host memory still works everywhere and is
    // staged through the same buffers; this only removes that copy.
    double *acquire_input();

    // --- accumulation ------------------------------------------------------

    // A += B, where B is column-major on the host: column j begins at
    // b + j*col_stride and its rows are contiguous (0 means tightly packed).
    // Throws std::domain_error on inf/NaN, and std::runtime_error if a column
    // was reserved too narrow or at too high an exponent for these values.
    //
    // **B must be page-locked and device-mapped** -- cudaHostAlloc with
    // cudaHostAllocMapped, cudaHostRegister with cudaHostRegisterMapped, or a
    // buffer from acquire_input(). Pageable memory is rejected, not copied.
    //
    // The library used to stage a copy of whatever it was handed, which cost
    // 1331 us at 65536x64 and made the fast path something a caller had to
    // know to ask for. Requiring pinned input removes the copy and, with it,
    // the two-contract problem: there is now one rule, the same for every
    // caller, and it is in the signature.
    //
    // That rule is the returned handle. The kernel streams B over PCIe and is
    // still reading it when this returns, so do not write to B until the
    // handle says the read has finished. Ignoring the handle is a statement
    // that B will not be written again.
    // `surveys` is optional, one entry per column, describing the same slice
    // this call reads. Supplied, it replaces the host survey this path would
    // otherwise take of B.
    void add_matrix_col_major(const double *b,
                                   const Survey *surveys = nullptr,
                                   std::size_t col_stride = 0);

    // A += B[0] + ... + B[count-1], all already resident in device memory, in
    // a single pass over the accumulation matrix.
    //
    // Identical results to submitting them one at a time. What changes is
    // traffic, and on this target that is the whole game: the accumulate is
    // bandwidth-bound at 97-99% of the card's streaming rate, moving ~8 bytes
    // of input and ~16*nlimbs of accumulation matrix per element. Batching turns
    // that into 8*count + 16*nlimbs, because every matrix's addend lands in the
    // same limbs of the same row and L1 absorbs the repeats.
    //
    // Blocked over the input set rather than held in registers on purpose. A
    // register array would have to be indexed at compile time, so the kernel
    // would be templated on the limb count -- which columns of one matrix do
    // not share, so one launch could not serve them.
    //
    // At most 16 matrices a call; they ride in the kernel's parameter block.
    //
    // Returns the same completion handle as everything else here, so a
    // producer can keep the buffers it is about to overwrite straight.
    // `surveys` is optional and shaped like `b`: surveys[k][j] for matrix k's
    // column j. Supplied, the survey kernel is not launched and its verdict is
    // not waited for, so this call stops synchronizing -- the 19.9 us drain
    // that otherwise sits between the submission and the accumulate.
    void add_matrices_col_major_device(const double *const *b,
                                            std::size_t count,
                                            const Survey *const *surveys
                                                = nullptr,
                                            std::size_t col_stride = 0);

    // A += B, for B already resident in device memory.
    //
    // The returned handle is what lets a producer cycle device buffers: fill
    // one while the accumulation matrix reads another, wait on the handle for
    // the one about to be overwritten, and never synchronize the stream.
    // Without it a caller's only recourse is synchronize(), which drains
    // everything and gives up exactly the overlap a rotation exists for.
    //
    // **Order your fill yourself.** This accumulation matrix's stream is
    // non-blocking, so nothing queued on the legacy null stream is implicitly
    // ordered against it -- including a cudaMemcpy from pageable host memory,
    // which returns before its DMA has landed. Two ways to be sure:
    //
    // Queue the fill on stream(), where the stream orders it against the
    // accumulate in both directions, and neither the handle nor a second
    // buffer is needed.
    //
    // Or use a stream of your own, synchronize it before submitting (that waits
    // for your copy, not for the accumulate), and wait on the returned handle
    // before refilling that buffer. Measured at 65536x64, a two-buffer rotation
    // costs 1301 us a batch that way against 1256 for the fill alone.
    // `surveys` is optional, one entry per column. Supplied, the survey
    // kernel is not launched and its verdict is not waited for, so this call
    // stops synchronizing.
    void add_matrix_col_major_device(const double *b,
                                          const Survey *surveys = nullptr,
                                          std::size_t col_stride = 0);

    // Blocks until every submitted accumulation has finished. Accumulation is
    // asynchronous: a call returns once the work is queued, so a caller
    // streaming batches keeps the device busy without doing anything special.
    // Readback synchronizes implicitly, so this is only needed for timing or
    // before reusing a caller-owned device buffer.
    void synchronize() const;

    // The stream this accumulation matrix launches on, for a producer that
    // would rather queue its own work there than coordinate two streams.
    // `CUstream_st *` is `cudaStream_t` spelled without the toolkit header,
    // so the result goes straight to any CUDA call that takes one.
    //
    // A submission is a kernel launch on this stream, so anything queued here
    // is ordered against it. That settles both halves of the buffer problem at
    // once: a fill queued before a submission is complete before the kernel
    // reads it, and a refill queued after one cannot begin until that kernel
    // retires. Neither needs a handshake, and neither touches the legacy null
    // stream, so a producer working this way never meets the serialization a
    // plain cudaMemcpy runs into.
    //
    // It is also how a caller learns when the device is done with a buffer it
    // submitted, which matters because the kernel is still reading that buffer
    // when the call returns -- a producer refilling it immediately corrupts
    // about 63% of entries, measured. Record an event here after submitting and
    // the ordering is the stream's:
    //
    //     acc.add_matrix_col_major_device(buf);
    //     cudaEventRecord(done, acc.stream());
    //     ...
    //     cudaEventSynchronize(done);   // or cudaEventQuery, to poll
    //
    // The library used to hand back a move-only handle that owned such an
    // event. It was removed once this accessor existed: the handle wrapped
    // three CUDA calls a caller can make directly, and made every submission
    // pay for an event whether or not anyone wanted one.
    //
    // The stream belongs to the accumulation matrix and is destroyed with it,
    // so do not destroy it, and remember that work queued here delays the
    // accumulates behind it as surely as they delay it.
    CUstream_st *stream() const;

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
    // device. The counterpart of AccumulationMatrix::entry_limbs.
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

        // The same derived width bound AccumulationMatrix keeps, and for the
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
    // costs, so the accumulation matrix work hides entirely behind the bus.
    // That only holds because the input is read once: the host survey is what
    // means the device never needs a second look at it.
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
    void locate(std::size_t &col, std::size_t &slot, std::size_t i,
                std::size_t j) const;
    void init_columns();
    void accumulate_device(const double *b, std::size_t col_stride,
                           const Survey *supplied);
    void launch_accumulate(const double *const *b, std::size_t count,
                           std::size_t col_stride);
    // Sentinels onto the device, then the verdict back. Every device-side
    // survey is bracketed by these two: launch begin_survey(), launch as many
    // survey kernels as there are matrices, then read what end_survey()
    // returns. The pointer it hands back is this object's own pinned buffer
    // and stays valid until the next survey.
    void begin_survey();
    const Survey *end_survey();
    void survey_device_inputs(const double *const *b, std::size_t count,
                              std::size_t col_stride,
                              const Survey *const *supplied);
    void validate_host_survey(const double *b, std::size_t col_stride,
                              const Survey *supplied);
    void require_all_reserved() const;
    void sync_descriptors();
    void harvest_occupancy();
    void reserve_from_extents(const int *low, const int *high,
                              const char *any, std::size_t count);
    // `count` is how many addends the extents cover, which is one per
    // matrix. A batch surveys its whole input set into a single pair of extents
    // and then fits the column once, so it has to say how many matrices that
    // pair stands for; the headroom term grows with the number of addends and
    // not with the range they span.
    void require_fit(std::size_t j, int min_exponent, int max_top,
                     std::size_t count = 1);
    void grow_column(std::size_t j, std::size_t needed);
    void rescale_column(std::size_t j, int new_exponent);
    void ensure_slots(std::size_t words);
    // Zeroes every limb array reserved since the last flush, in one launch.
    // Const because it completes deferred work rather than changing what the
    // accumulation matrix holds -- the limbs read as zero either way, and
    // readback paths are const and must be able to complete it.
    void flush_pending_zero() const;

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
    bool owns_stream_ = true;  // false when the caller supplied it
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

    // The same pair for the detector's verdict, brought across by the same
    // elected block in the same pass.
    unsigned *flags_device_ = nullptr;
    unsigned *flags_host_ = nullptr;

    void *descriptors_ = nullptr;  // device ColumnDesc[]
    bool descriptors_stale_ = true;
    // Page-locked staging for the descriptor upload, and an event marking when
    // the device has finished reading it. cudaMemcpyAsync from *pageable*
    // memory synchronizes the stream before it starts, so copying out of an
    // ordinary vector would drain the pipeline however asynchronous it looks.
    // Limb arrays allocated but not yet zeroed. reserve_column issues one
    // allocation per limb position, and zeroing each with its own
    // cudaMemsetAsync costs far more than one kernel over all of them: at 512
    // arrays, 827 us of memsets against 59 for a single launch. Deferred to
    // the first launch or readback so a whole pre-sized accumulation matrix is
    // zeroed at once, which is where the count is largest.
    mutable std::vector<limbs::limb_t *> zero_ptr_;
    mutable std::vector<std::size_t> zero_rows_;
    mutable void *zero_targets_ = nullptr;  // device ZeroTarget[]
    mutable std::size_t zero_capacity_ = 0;

    void *desc_host_ = nullptr;  // pinned ColumnDesc[]
    CUevent_st *ev_desc_ = nullptr;
    bool desc_in_flight_ = false;
    void *survey_out_ =
        nullptr;  // device survey results, for device-side input
    // Where those results are read back to. Pinned, so the copy out is a DMA
    // rather than one the driver stages through a buffer of its own; not
    // mapped, because the kernel reduces into `survey_out_` with atomics and
    // those belong in device memory.
    void *survey_host_ = nullptr;  // pinned Survey[]

    // Per-column stored length and first row, device-side. Set at
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

}  // namespace truesum
