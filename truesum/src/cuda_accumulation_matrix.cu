#include <cuda_runtime.h>

#include <algorithm>
#include <climits>
#include <cstring>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "truesum/accumulation_matrix.hpp"
#include "truesum/cuda_accumulation_matrix.hpp"
#include "kernels.hpp"

namespace truesum {
namespace {

using limbs::limb_t;

// The header forward-declares these rather than including cuda_runtime.h. If
// CUDA ever renamed the underlying structs this would stop compiling here,
// which is the point.
static_assert(std::is_same<cudaStream_t, CUstream_st *>::value,
              "cudaStream_t is no longer CUstream_st*");
static_assert(std::is_same<cudaEvent_t, CUevent_st *>::value,
              "cudaEvent_t is no longer CUevent_st*");

[[noreturn]] void
cuda_fail(cudaError_t status, const char *call, const char *file, int line)
{
    std::ostringstream os;
    os << "truesum: " << call << " failed at " << file << ":" << line << ": "
       << cudaGetErrorName(status) << " -- " << cudaGetErrorString(status);
    throw std::runtime_error(os.str());
}

inline void
cuda_check(cudaError_t status, const char *call, const char *file, int line)
{
    if (cudaSuccess != status) cuda_fail(status, call, file, line);
}

// Every checked runtime call reads exactly like the call it makes, with one
// extra parenthesis: cuda(Malloc(&p, bytes)) is cudaMalloc(&p, bytes). The
// macro pastes the prefix back on for the call and stringizes it for the
// diagnostic, so a failure reports the function rather than an opaque
// expression. A blocking wait is cuda(StreamSynchronize(s)), which pastes the
// same way, so there is one mechanism rather than several.
//
// Launches are not separately checked. Every launch parameter here is either a
// compile-time constant or bounded at construction, so cudaGetLastError after
// one has nothing to report; a genuine fault is sticky and surfaces at the
// next stream synchronize, which readback does implicitly.
//
// Being function-like, the macro only fires on `cuda` followed by `(`, which
// leaves every ordinary cudaXxx token -- cudaError_t, cudaSuccess, the
// unchecked cudaFree below -- alone. It is #undef'd at the end of the file.
//
// Failures throw rather than jumping to a cleanup label: this is C++, the
// device allocations are owned by the object, and the rest of the library
// already reports misuse by throwing. The destructor is the only place that
// swallows a status, because it must not throw.
#define cuda(call) cuda_check(cuda##call, "cuda" #call, __FILE__, __LINE__)

constexpr unsigned kBlock = 256;  // a power of two; the reductions rely on it

// Mirrors kernels.hpp, which this file cannot include from device code.
constexpr unsigned kBadNonFinite = 1;
constexpr unsigned kBadExponent = 2;
constexpr unsigned kBadWidth = 4;
constexpr unsigned kBadOverflow = 8;

// Blocks in the whole grid, not just its x extent. Both kernels grid-stride,
// so any smaller grid stays correct -- it only gives each thread more rows.
//
// Capping the total is what makes the survey's reduction pay for itself. Left
// uncapped, 65536 rows over 64 columns launches 16384 blocks, each reducing a
// single value per thread and then paying a full eight-step shared-memory tree
// to do it. Capped, each thread reduces sixteen rows serially first and the tree
// is amortized across them: measured 195 -> 113 us at that shape, and
// 746 -> 423 us at 262144 rows.
constexpr unsigned kMaxBlocks = 1024;

// A grid that covers `rows` a thread at a time, then shrunk to kMaxBlocks.
// Subsumes the old separate cap on the x extent: with one column the two are
// the same bound.
dim3
launch_grid(std::size_t rows, std::size_t cols)
{
    unsigned gx = static_cast<unsigned>((rows + kBlock - 1) / kBlock);
    const unsigned per_column =
        0 == cols ? kMaxBlocks : static_cast<unsigned>(kMaxBlocks / cols);
    if (gx > per_column) gx = per_column;
    if (0 == gx) gx = 1;
    return dim3(gx, static_cast<unsigned>(cols));
}

// The survey kernel reduces straight into the public Survey. Its two extents
// are ints at offsets 0 and 4, and the struct's 12-byte size keeps every entry
// 4-aligned, so both take atomics directly. The bools need none: every writer
// writes true, so the race is benign and a plain store is correct.
//
// There is deliberately no device-side twin of this struct. One existed, with
// int in place of the bools, and it bought nothing but a conversion at every
// boundary.

// Sentinels for the reduction, and the tidy-up afterwards. Both are one thread
// a column and cost a launch each; they exist so a caller never sees a
// sentinel and cannot tell which side computed the survey.
__global__ void
survey_init_kernel(Survey *__restrict__ out, unsigned ncols)
{
    const unsigned j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= ncols) return;
    out[j].min_exponent = INT_MAX;
    out[j].max_top = INT_MIN;
    out[j].any = false;
    out[j].nonfinite = false;
}

__global__ void
survey_finish_kernel(Survey *__restrict__ out, unsigned ncols)
{
    const unsigned j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= ncols) return;
    // An empty or non-finite column reports zero extents, as the host survey
    // does -- the sentinels above are not a value any caller should see.
    if (!out[j].any || out[j].nonfinite) {
        out[j].min_exponent = 0;
        out[j].max_top = 0;
    }
}

// A column's stored length and where it starts. Set at construction, so
// unlike ColumnDesc this never goes stale and both kernels can read it
// without being sequenced against a descriptor sync.
struct ColumnShape {
    unsigned rows;
    unsigned first_row;
};

// The matrices one launch batches together. Passed by value, so it rides in the
// kernel parameter block (constant memory) and needs no allocation or copy.
constexpr unsigned kMaxBatchInputs = 16;
struct InputSet {
    const double *p[kMaxBatchInputs];
    unsigned n;
};

struct ColumnDesc {
    // One base per limb position, so widening appends rather than moves.
    // A warp's threads all work different rows of the same limb column, so
    // this pointer load is a broadcast of one value, L1-resident for the life
    // of the kernel: the array is nlimbs * 8 bytes, 128 for a 16-limb column.
    limb_t *const *bases;
    int exponent;
    unsigned nlimbs;
};

// Every kernel below walks its rows with a grid-stride loop, spelled out at
// each one so the thread and block indices stay visible -- they are what the
// loop is about. The widening cast is not decoration: blockIdx.x, blockDim.x
// and gridDim.x are all unsigned int, and both products overflow 32 bits above
// 4.2M threads. launch_grid caps well below that today, which is exactly why
// the cast is easy to drop and worth keeping.

// The device twin of the scalar kernel's `split`, field for field. Kept in
// lockstep with it deliberately: the two implementations have to agree about
// what a double means before they can be compared bit for bit.
//
//   normal    (biased != 0): m = 2^52 | frac, e = biased - 1075
//   denormal (biased == 0): m = frac,        e = -1074
//
// Both collapse to e = max(biased, 1) - 1075.
__device__ inline void
split_device(double v, unsigned long long &mantissa, int &exponent, int &top,
             bool &negative, bool &nonfinite)
{
    const unsigned long long bits =
        static_cast<unsigned long long>(__double_as_longlong(v));
    const unsigned long long biased = (bits >> 52) & 0x7FFull;
    const unsigned long long frac = bits & ((1ull << 52) - 1);

    unsigned long long m = frac;
    if (0 != biased) m |= 1ull << 52;
    int e = static_cast<int>(biased < 1 ? 1 : biased) - 1075;

    negative = (bits >> 63) != 0;
    nonfinite = (biased == 0x7FFull);
    // Normalizing to odd raises the exponent by the trailing zero count and
    // lowers the significand's width by the same amount, so `top` is
    // unaffected and never needs that count.
    top = e + (0 == m ? 0 : 64 - __clzll(static_cast<long long>(m)));
    if (0 != m) {
        const int tz = __ffsll(static_cast<long long>(m)) - 1;
        m >>= tz;
        e += tz;
    }
    mantissa = m;
    exponent = e;
}

// The CPU survey splits into two passes so the significand can be skipped once
// a column's scale is set. That does not pay here: this pass is bound by
// reading the column, and the trailing-zero count is a few ALU ops on a value
// already in registers. One exact pass is both simpler and cheaper.
__global__ void
survey_kernel(const double *__restrict__ values,
              const ColumnShape *__restrict__ shapes, unsigned uniform_rows,
              std::size_t col_stride, Survey *__restrict__ out)
{
    // Only the extents need a tree. `any` is not reduced at all -- it is
    // implied by whether the minimum is still its sentinel -- and `nonfinite`
    // is written directly by whichever threads see one, which in the ordinary
    // case is none of them. Two arrays instead of four, so 2 KB a block rather
    // than 4, and half the work in every step of the tree below.
    __shared__ int s_min[kBlock];
    __shared__ int s_max[kBlock];

    const unsigned tid = threadIdx.x;
    const unsigned j = blockIdx.y;
    // A null `shapes` means every column has the same shape: `uniform_rows`
    // tall and starting at row 0, which is what surveying a whole matrix looks
    // like as opposed to an accumulation matrix's stored triangle. The shape
    // then rides in the parameter block, so the standalone survey below needs
    // no device allocation for a table of identical entries.
    const ColumnShape sh = shapes ? shapes[j] : ColumnShape{uniform_rows, 0};
    const double *col = values + static_cast<std::size_t>(j) * col_stride +
                        sh.first_row;
    const std::size_t rows = sh.rows;

    int tmin = INT_MAX;
    int tmax = INT_MIN;
    int tbad = 0;

    for (std::size_t i = (std::size_t)blockIdx.x * blockDim.x + threadIdx.x;
                     i < rows;
                     i += (std::size_t)blockDim.x * gridDim.x) {
        unsigned long long m;
        int e, top;
        bool neg, bad;
        split_device(col[i], m, e, top, neg, bad);
        if (bad) tbad = 1;
        if (0 != m) {
            if (e < tmin) tmin = e;
            if (top > tmax) tmax = top;
        }
    }

    s_min[tid] = tmin;
    s_max[tid] = tmax;
    // Rare, and unreduced: a plain store of true from every thread that saw a
    // non-finite value. They all write the same value, so the race is benign,
    // and a clean batch performs no store at all.
    if (tbad) out[j].nonfinite = true;
    __syncthreads();
    for (unsigned s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            if (s_min[tid + s] < s_min[tid]) s_min[tid] = s_min[tid + s];
            if (s_max[tid + s] > s_max[tid]) s_max[tid] = s_max[tid + s];
        }
        __syncthreads();
    }

    // One atomic per block rather than one per thread. A minimum still at its
    // sentinel means no lane in this block was live, which is exactly what an
    // `any` reduction would have told us -- no biased exponent can reach
    // INT_MAX, so the sentinel is unambiguous.
    if (0 == tid && INT_MAX != s_min[0]) {
        atomicMin(&out[j].min_exponent, s_min[0]);
        atomicMax(&out[j].max_top, s_max[0]);
        out[j].any = true;
    }
}

// Thread i owns row i, so no two threads touch the same limb word and the
// accumulation needs no atomics. Consecutive threads read consecutive words of
// limbs[p], which is what limb-major buys on this target: a warp's 32 accesses
// are 256 contiguous bytes, fully coalesced.
//
// Unlike the AVX-512 kernel, each thread walks its *own* limb range starting
// at its own offset -- there is no shared loop counter, so a warp's trip count
// is the maximum over its threads rather than the range of their offsets.
// Divergent exponents cost coalescing here (threads land in different limb
// arrays in the same instruction), not instruction count.
template <int kRadix>
__global__ void
accumulate_kernel(const ColumnDesc *__restrict__ cols, InputSet in,
                  const ColumnShape *__restrict__ shapes,
                  std::size_t col_stride, int *__restrict__ occupancy_device,
                  int *__restrict__ occupancy_host,
                  unsigned *__restrict__ flags_device,
                  unsigned *__restrict__ flags_host,
                  unsigned *__restrict__ ticket, unsigned ncols)
{
    // The addend split, the offset arithmetic and the limb mask below are all
    // written in terms of kRadix, so a reduced radix would drop in here. It is
    // not worth doing on this target: deleting the carry chain outright --
    // the most carry-save could ever save -- measured -2.4% to +2.8% on the
    // device, because this kernel is bandwidth-bound and the arithmetic hides
    // behind the memory. The same experiment on AVX-512 is worth 25-52%. An
    // earlier version of this comment called carry-save the device's next
    // milestone, which pointed at the one target where it was measured not to
    // pay. See docs/simd-design.md.
    static_assert(kRadix == 64,
                  "only the canonical radix is implemented");
    constexpr unsigned long long kLimbMask =
        kRadix == 64 ? ~0ull : ((1ull << kRadix) - 1);

    const unsigned j = blockIdx.y;
    const ColumnDesc c = cols[j];
    const ColumnShape sh = shapes[j];
    const std::size_t col_off =
        static_cast<std::size_t>(j) * col_stride + sh.first_row;
    const std::size_t rows = sh.rows;

    // Only the occupancy needs a tree. The flags are ORed straight into the
    // per-column staging by whichever threads raise one, which in the ordinary
    // case is none: a clean batch performs no atomic and needs no array.
    __shared__ int s_lim[kBlock];
    const unsigned tid = threadIdx.x;
    // Highest limb position this thread disturbs. The carry loop already knows
    // where it stopped, so this costs a comparison rather than a pass.
    int t_lim = -1;
    // Where this batch contradicts what the column was sized for. Pre-sized
    // from producer metadata there is no survey to catch it, and the loop
    // below would otherwise drop such a value without a word: a negative
    // shift converts to an enormous `off` and the carry loop simply never
    // runs. Recorded, not acted on -- the arithmetic is left exactly as it
    // was so the checks cost three comparisons and no divergence.
    unsigned t_flags = 0;
    unsigned long long t_overflow = 0;

    for (std::size_t i = (std::size_t)blockIdx.x * blockDim.x + threadIdx.x;
                     i < rows;
                     i += (std::size_t)blockDim.x * gridDim.x) {
        // Blocking over the input set, not over the accumulation matrix. Every
        // matrix's addend lands in the same limbs of the same row, back to
        // back, so those read-modify-writes hit L1 and only the first read and
        // the last write reach DRAM. Traffic an element goes from 8 + 16*nlimbs
        // to 8*n + 16*nlimbs, and it needs no register array -- so the limb
        // count does not have to be a compile-time constant, which it could not
        // be: columns of one matrix have their own widths.
        for (unsigned bi = 0; bi < in.n; ++bi) {
            const double *col = in.p[bi] + col_off;
            unsigned long long m;
            int e, top;
            bool neg, bad;
            split_device(col[i], m, e, top, neg, bad);
            if (bad) t_flags |= kBadNonFinite;
            if (0 == m) continue;

            // A negative shift converts to an `off` far above nlimbs, so the
            // loop below does not run and the value is silently dropped. That
            // is the failure the flags exist to report.
            const int shift = e - c.exponent;
            const unsigned off = static_cast<unsigned>(shift) / kRadix;
            if (shift < 0) t_flags |= kBadExponent;
            // The addend's top must stay below the sign bit, which also keeps
            // the sign test at the top limb exact.
            if (top - c.exponent >= static_cast<int>(kRadix * c.nlimbs)) {
                t_flags |= kBadWidth;
            }
            const unsigned bit = static_cast<unsigned>(shift) % kRadix;

            // A 53-bit significand at intra-limb offset `bit` spans two limbs
            // at radix 64 and at radix 52 alike.
            unsigned long long lo, hi;
            if (0 == bit) {
                lo = m & kLimbMask;
                hi = kRadix == 64 ? 0ull : (m >> kRadix) & kLimbMask;
            } else {
                lo = (m << bit) & kLimbMask;
                hi = (m >> (kRadix - bit)) & kLimbMask;
            }

            // Add and subtract share the loop shape; `carry` is a borrow when
            // the lane is negative. Past the addend with nothing propagating,
            // no higher limb can change.
            unsigned long long carry = 0;
            for (unsigned p = off; p < c.nlimbs; ++p) {
                const unsigned long long a =
                    (p == off) ? lo : ((p == off + 1) ? hi : 0ull);
                limb_t *dst = c.bases[p] + i;
                const unsigned long long x = *dst;
                unsigned long long written;
                // Signed overflow at the top limb -- a negative value that lost
                // its sign bit, or a non-negative one that gained it -- kept as
                // a bit rather than branched on, so lanes do not diverge on the
                // sign of a running sum.
                if (neg) {
                    const unsigned long long d = x - a;
                    const unsigned long long b1 = (x < a) ? 1ull : 0ull;
                    const unsigned long long d2 = d - carry;
                    const unsigned long long b2 = (d < carry) ? 1ull : 0ull;
                    written = d2;
                    carry = b1 | b2;
                    if (p + 1 == c.nlimbs) t_overflow |= (x & ~d2) >> 63;
                } else {
                    const unsigned long long s = x + a;
                    const unsigned long long c1 = (s < x) ? 1ull : 0ull;
                    const unsigned long long s2 = s + carry;
                    const unsigned long long c2 = (s2 < s) ? 1ull : 0ull;
                    written = s2;
                    carry = c1 | c2;
                    if (p + 1 == c.nlimbs) t_overflow |= (~x & s2) >> 63;
                }
                *dst = written;

                // Significant, not merely written. All-zeros and all-ones are
                // exactly what a two's complement sign extension leaves behind,
                // and a borrow out of a negative addend writes all-ones every
                // limb to the top of the column -- so counting writes would
                // report the full width for any column that ever goes negative.
                // A value's topmost significant limb can never be all-ones when
                // positive (the sign bit would be set) nor all-zeros when
                // negative, so this test finds exactly that limb.
                if (0ull != written && ~0ull != written &&
                    static_cast<int>(p) > t_lim) {
                    t_lim = static_cast<int>(p);
                }
                if (0 == carry && p >= off + 1) break;
            }
      }
    }

    // Batch the per-thread maxima, reduce across blocks in device memory with
    // ordinary atomics, then elect one block to copy the finished array to
    // host memory. Atomics never cross PCIe: measured, that costs 260x.
    if (0 != t_overflow) t_flags |= kBadOverflow;
    if (0 != t_flags) atomicOr(&flags_device[j], t_flags);

    s_lim[tid] = t_lim;
    __syncthreads();
    for (unsigned s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            if (s_lim[tid + s] > s_lim[tid]) s_lim[tid] = s_lim[tid + s];
        }
        __syncthreads();
    }
    if (0 == tid && s_lim[0] >= 0) atomicMax(&occupancy_device[j], s_lim[0]);

    // Every block must publish before the elected one reads.
    __threadfence();
    __shared__ bool last;
    if (0 == tid) {
        last = (atomicAdd(ticket, 1u) == gridDim.x * gridDim.y - 1);
    }
    __syncthreads();
    if (!last) return;

    // One block, one thread per column: the write that does cross PCIe, and
    // the only one.
    //
    // The staging is deliberately *not* re-armed. Occupancy is a high-water
    // mark over every batch, so leaving it to accumulate on the device is both
    // what the fit test wants and what keeps the host from having to observe
    // each launch -- resetting it here meant a reader saw only the last
    // batch's maximum, which lost limbs that earlier batches had needed.
    for (unsigned k = tid; k < ncols; k += blockDim.x) {
        occupancy_host[k] = occupancy_device[k];
        flags_host[k] = flags_device[k];
    }
    if (0 == tid) *ticket = 0;
}

// One limb array to be zeroed, and the kernel that does the lot.
struct ZeroTarget {
    limb_t *p;
    unsigned rows;
};

__global__ void
zero_kernel(const ZeroTarget *__restrict__ targets)
{
    const ZeroTarget z = targets[blockIdx.y];
    for (std::size_t i = (std::size_t)blockIdx.x * blockDim.x + threadIdx.x;
                     i < z.rows;
                     i += (std::size_t)blockDim.x * gridDim.x) {
        z.p[i] = 0;
    }
}

// A newly appended limb array holds the sign extension of the one below it:
// all-ones where that limb is negative, all-zeros otherwise. This is the whole
// cost of widening in a limb-major layout -- nothing already allocated moves.
__global__ void
sign_fill_kernel(limb_t *__restrict__ dst, const limb_t *__restrict__ src,
                 std::size_t rows)
{
    for (std::size_t i = (std::size_t)blockIdx.x * blockDim.x + threadIdx.x;
                     i < rows;
                     i += (std::size_t)blockDim.x * gridDim.x) {
        dst[i] = static_cast<limb_t>(static_cast<long long>(src[i]) >> 63);
    }
}

// Shifts every entry of a column left by `shift` bits, in place.
//
// In place is safe because each thread owns one row, so the ordering that
// matters is within a thread rather than across them. Walking limb positions
// downward, position k reads k-word and k-word-1, both at or below k: strictly
// below when word >= 1, and when word == 0 the read of k happens before its
// own write. Nothing a thread still needs has been overwritten.
//
// The caller widens first, so every position in [0, nlimbs) already holds a
// valid sign-extended limb and there is no special case above the old top.
__global__ void
shift_left_kernel(limb_t *const *__restrict__ bases, unsigned nlimbs,
                  std::size_t rows, unsigned shift)
{
    const unsigned word = shift / 64;
    const unsigned bit = shift % 64;
    for (std::size_t i = (std::size_t)blockIdx.x * blockDim.x + threadIdx.x;
                     i < rows;
                     i += (std::size_t)blockDim.x * gridDim.x) {
        for (int k = static_cast<int>(nlimbs) - 1; k >= 0; --k) {
            const int hi = k - static_cast<int>(word);
            limb_t v = 0;
            if (hi >= 0) {
                v = bases[hi][i];
                if (0 != bit) {
                    v <<= bit;
                    if (hi >= 1) v |= bases[hi - 1][i] >> (64 - bit);
                }
            }
            bases[k][i] = v;
        }
    }
}

std::size_t
ceil_log2(std::size_t n)
{
    std::size_t b = 0;
    while ((std::size_t{1} << b) < n) ++b;
    return b;
}

std::size_t
limbs_for_bits(std::size_t bits)
{
    return (bits + limbs::kLimbBits - 1) / limbs::kLimbBits;
}

// Both pointers a caller hands the standalone survey are dereferenced by the
// kernel and by nothing on the host, so each has to be memory the device can
// reach. Pageable memory faults inside the launch, away from the mistake;
// one driver query on a setup path buys a diagnostic that says what to
// allocate instead. What comes back is the device alias, because a mapped host
// allocation is not obliged to share its address with the device.
void *
device_alias(const void *p, const char *what)
{
    cudaPointerAttributes attr{};
    const cudaError_t st = cudaPointerGetAttributes(&attr, p);
    if (cudaSuccess != st || nullptr == attr.devicePointer) {
        cudaGetLastError();  // an unregistered pointer leaves this sticky
        throw std::invalid_argument(
            std::string("truesum: ") + what +
            " from the device, so it must be device memory or page-locked "
            "mapped host memory -- cudaMalloc, cudaHostAlloc with "
            "cudaHostAllocMapped, or cudaHostRegister with "
            "cudaHostRegisterMapped");
    }
    return attr.devicePointer;
}

}  // namespace

void
survey_matrix_col_major_device(Survey *out, const double *b,
                               std::size_t rows, std::size_t cols,
                               std::size_t col_stride, CUstream_st *stream)
{
    if (0 == rows || 0 == cols) return;
    const std::size_t stride = col_stride ? col_stride : rows;

    // The kernel writes `out` and reads `b`, and the host touches neither, so
    // both are device pointers whatever kind of memory the caller allocated.
    Survey *dst = static_cast<Survey *>(
        device_alias(out, "survey_matrix_col_major_device writes `out`"));
    const double *src = static_cast<const double *>(
        device_alias(b, "survey_matrix_col_major_device reads `b`"));

    // Every column is the full height, since this surveys a matrix rather than
    // an accumulation matrix's stored triangle. That uniform shape goes in the
    // parameter block, so the survey allocates no device memory at all: a
    // caller who owns both buffers is entitled to a call that neither allocates
    // behind their back nor perturbs the allocator between their own launches.
    const unsigned n = static_cast<unsigned>(cols);
    const unsigned init_blocks = (n + kBlock - 1) / kBlock;
    survey_init_kernel<<<init_blocks, kBlock, 0, stream>>>(dst, n);
    survey_kernel<<<launch_grid(rows, cols), kBlock, 0, stream>>>(
        src, nullptr, static_cast<unsigned>(rows), stride, dst);
    survey_finish_kernel<<<init_blocks, kBlock, 0, stream>>>(dst, n);
    // No wait. The three launches are ordered against each other by the stream
    // they share, and when the survey is readable is the caller's business:
    // they may want it on the device, or beside their own launches on their own
    // stream, and either way they can synchronize or record an event.
}

bool
cuda_available()
{
    int n = 0;
    return cudaSuccess == cudaGetDeviceCount(&n) && n > 0;
}

CudaAccumulationMatrix::CudaAccumulationMatrix(std::size_t rows,
                                               std::size_t cols,
                                               CUstream_st *stream)
    : rows_(rows), cols_(cols), cols_state_(cols)
{
    if (!cuda_available()) {
        throw std::runtime_error("truesum: no usable CUDA device");
    }
    // Mapping has to be enabled before the context exists, so this fails
    // harmlessly if one is already active with the flag set. The cudaHostAlloc
    // below is the check that matters -- it fails loudly if mapping is really
    // unavailable.
    cudaSetDeviceFlags(cudaDeviceMapHost);
    // Columns become the grid's y dimension, the only launch parameter here
    // that is not constant at compile time or clamped. Bounding it once, up front,
    // is what lets every launch below go unchecked.
    int max_grid_y = 0;
    cuda(DeviceGetAttribute(&max_grid_y, cudaDevAttrMaxGridDimY, 0));
    if (cols_ > static_cast<std::size_t>(max_grid_y)) {
        std::ostringstream os;
        os << "truesum: " << cols_ << " columns exceeds this device's grid y "
           << "limit of " << max_grid_y;
        throw std::runtime_error(os.str());
    }
    // The caller's stream if there is one, and it stays theirs: not destroyed
    // here, and carrying whatever flags they gave it.
    //
    // Ours is blocking, which cudaStreamCreate gives by default and which is
    // the safe choice rather than the fast one. A blocking stream is implicitly
    // ordered against the legacy null stream, so a caller who fills a device
    // buffer with cudaMemcpy and submits it cannot be caught out -- and that is
    // worth protecting, because a cudaMemcpy from *pageable* host memory
    // returns once the source is staged, with the DMA still in flight. Made
    // non-blocking, this suite failed 3 runs in 40 on exactly that shape.
    //
    // The cost is the overlap: a producer sharing the null stream waits for the
    // accumulate, 1776 us a batch against 1301 at 65536x64. A caller who wants
    // that back supplies a stream of their own and takes on the ordering.
    if (nullptr != stream) {
        st_compute_ = stream;
        owns_stream_ = false;
    } else {
        cuda(StreamCreate(&st_compute_));
    }

    if (0 != cols_) {
        cuda(Malloc(&occupancy_device_, cols_ * sizeof(int)));
        cuda(Memset(occupancy_device_, 0xFF, cols_ * sizeof(int)));  // -1
        cuda(HostAlloc(&occupancy_host_, cols_ * sizeof(int),
                       cudaHostAllocMapped));
        for (std::size_t j = 0; j < cols_; ++j) occupancy_host_[j] = -1;
        cuda(Malloc(&flags_device_, cols_ * sizeof(unsigned)));
        cuda(Memset(flags_device_, 0, cols_ * sizeof(unsigned)));
        cuda(HostAlloc(&flags_host_, cols_ * sizeof(unsigned),
                       cudaHostAllocMapped));
        for (std::size_t j = 0; j < cols_; ++j) flags_host_[j] = 0;
        cuda(Malloc(&ticket_, sizeof(unsigned)));
        cuda(Memset(ticket_, 0, sizeof(unsigned)));
        cuda(Malloc(&descriptors_, cols_ * sizeof(ColumnDesc)));
        cuda(HostAlloc(&desc_host_, cols_ * sizeof(ColumnDesc),
                       cudaHostAllocDefault));
        cuda(EventCreateWithFlags(&ev_desc_, cudaEventDisableTiming));
        cuda(Malloc(&survey_out_, cols_ * sizeof(Survey)));
        cuda(HostAlloc(&survey_host_, cols_ * sizeof(Survey),
                       cudaHostAllocDefault));
    }
    init_columns();
}

CudaAccumulationMatrix::CudaAccumulationMatrix(std::size_t n, Uplo uplo,
                                               CUstream_st *stream)
    : CudaAccumulationMatrix(n, n, stream)
{
    symmetric_ = true;
    uplo_ = uplo;
    init_columns();
}

// Fills in each column's stored length and first row, and publishes them to
// the device. Called again by the symmetric constructor because the
// delegated-to one has already run with the full shape.
void
CudaAccumulationMatrix::init_columns()
{
    std::vector<ColumnShape> shapes(cols_);
    for (std::size_t j = 0; j < cols_; ++j) {
        Column &c = cols_state_[j];
        if (symmetric_) {
            c.rows = Uplo::Lower == uplo_ ? rows_ - j : j + 1;
            c.first_row = Uplo::Lower == uplo_ ? j : 0;
        } else {
            c.rows = rows_;
            c.first_row = 0;
        }
        shapes[j].rows = static_cast<unsigned>(c.rows);
        shapes[j].first_row = static_cast<unsigned>(c.first_row);
    }
    if (0 == cols_) return;
    if (nullptr == shapes_) {
        cuda(Malloc(&shapes_, cols_ * sizeof(ColumnShape)));
    }
    cuda(Memcpy(shapes_, shapes.data(), cols_ * sizeof(ColumnShape),
                cudaMemcpyHostToDevice));
}

CudaAccumulationMatrix::CudaAccumulationMatrix(
    std::size_t rows, std::size_t cols,
    const std::vector<std::vector<Survey>> &surveys, CUstream_st *stream)
    : CudaAccumulationMatrix(rows, cols, stream)
{
    reserve_from_surveys(surveys);
}

CudaAccumulationMatrix::CudaAccumulationMatrix(
    std::size_t n, Uplo uplo, const std::vector<std::vector<Survey>> &surveys,
    CUstream_st *stream)
    : CudaAccumulationMatrix(n, uplo, stream)
{
    reserve_from_surveys(surveys);
}

// The same aggregate AccumulationMatrix forms: minimum of the minima, maximum
// of the maxima, count from the outer size. Reducing to one reservation per
// column is what makes submission order irrelevant -- any matrix inside these
// extents fits, so the accumulation matrix never needs to know which one it is
// being handed.
//
// Unlike the CPU, a column with nothing in it is still reserved: the device
// can grow a column but only from the host between launches, and a pre-sized
// accumulation matrix is meant never to go back to the host at all.
void
CudaAccumulationMatrix::reserve_from_surveys(
    const std::vector<std::vector<Survey>> &surveys)
{
    if (surveys.empty()) {
        throw std::invalid_argument(
            "truesum: pre-sizing needs the surveys of at least one matrix");
    }
    for (const auto &one : surveys) {
        if (one.size() != cols_) {
            throw std::invalid_argument(
                "truesum: each survey must have one entry per column");
        }
    }

    const std::size_t headroom = ceil_log2(surveys.size() + 1) + 1;
    for (std::size_t j = 0; j < cols_; ++j) {
        bool any = false;
        int low = 0, high = 0;
        for (const auto &one : surveys) {
            if (!one[j].any) continue;
            if (!any) {
                low = one[j].min_exponent;
                high = one[j].max_top;
                any = true;
            } else {
                low = std::min(low, one[j].min_exponent);
                high = std::max(high, one[j].max_top);
            }
        }
        if (!any) {
            reserve_column(j, 0, limbs::kLimbBits);
            continue;
        }
        reserve_column(j, low,
                       static_cast<std::size_t>(high - low) + headroom);
    }

    presized_ = true;
    declared_matrices_ = surveys.size();
}
// The same sizing without the promise, and on this target it is also how a
// caller reserves columns from metadata rather than from data. Pre-sizing here
// removes the survey kernel and its round trip, which this does not: the
// accumulation matrix still surveys what it is handed. What it removes is the
// rescale, which on the device is not merely expensive but forbidden mid-launch,
// so a column that starts low enough is a column that cannot be caught out.
void
CudaAccumulationMatrix::reserve_for_surveys(
    const std::vector<std::vector<Survey>> &surveys)
{
    const bool was_presized = presized_;
    const std::size_t declared = declared_matrices_;
    reserve_from_surveys(surveys);
    presized_ = was_presized;
    declared_matrices_ = declared;
}


std::size_t
CudaAccumulationMatrix::column_rows(std::size_t j) const
{
    if (j >= cols_)
        throw std::out_of_range("truesum: column index out of range");
    return cols_state_[j].rows;
}

std::size_t
CudaAccumulationMatrix::column_first_row(std::size_t j) const
{
    if (j >= cols_)
        throw std::out_of_range("truesum: column index out of range");
    return cols_state_[j].first_row;
}

// The same batch AccumulationMatrix uses: Lower puts (i, j) and (j, i) both in
// column min(i, j) at slot |i - j|, Upper in column max(i, j) at slot
// min(i, j).
void
CudaAccumulationMatrix::locate(std::size_t &col, std::size_t &slot,
                               std::size_t i, std::size_t j) const
{
    if (!symmetric_) {
        col = j;
        slot = i;
        return;
    }
    if (Uplo::Lower == uplo_) {
        col = std::min(i, j);
        slot = i < j ? j - i : i - j;
    } else {
        col = std::max(i, j);
        slot = std::min(i, j);
    }
}

CudaAccumulationMatrix::~CudaAccumulationMatrix()
{
    // Deliberately unchecked: a destructor must not throw, and there is
    // nothing useful to do about a failed free during teardown.
    //
    // Drained first. The limb arrays are released stream-ordered, and a kernel
    // still reading them has to have retired before that release can recycle
    // the memory. The blocking default would have covered this; a caller's
    // stream does not have to be blocking.
    if (nullptr != st_compute_) (void)cudaStreamSynchronize(st_compute_);
    for (auto &c : cols_state_) {
        for (limb_t *p : c.bases) (void)cudaFreeAsync(p, st_compute_);
        cudaFree(c.dev_bases);
    }
    cudaFree(descriptors_);
    cudaFree(zero_targets_);
    cudaFreeHost(desc_host_);
    if (nullptr != ev_desc_) (void)cudaEventDestroy(ev_desc_);
    cudaFree(survey_out_);
    cudaFreeHost(survey_host_);
    cudaFree(shapes_);
    cudaFree(occupancy_device_);
    cudaFreeHost(occupancy_host_);
    cudaFree(flags_device_);
    cudaFreeHost(flags_host_);
    cudaFree(ticket_);
    for (Slot &s : slots_) {
        if (s.mapped) cudaFreeHost(s.mapped);
        if (s.ev_done) (void)cudaEventDestroy(s.ev_done);
    }
    if (st_compute_) {
        if (owns_stream_) (void)cudaStreamDestroy(st_compute_);
    }
}

void
CudaAccumulationMatrix::check_index(std::size_t i, std::size_t j) const
{
    if (i >= rows_ || j >= cols_) {
        throw std::out_of_range("truesum: matrix index out of range");
    }
}

int
CudaAccumulationMatrix::column_exponent(std::size_t j) const
{
    if (j >= cols_)
        throw std::out_of_range("truesum: column index out of range");
    return cols_state_[j].exponent;
}

std::size_t
CudaAccumulationMatrix::column_limbs(std::size_t j) const
{
    if (j >= cols_)
        throw std::out_of_range("truesum: column index out of range");
    return cols_state_[j].nlimbs;
}

void
CudaAccumulationMatrix::reserve_column(std::size_t j, int exponent,
                                       std::size_t bits)
{
    if (j >= cols_)
        throw std::out_of_range("truesum: column index out of range");
    Column &c = cols_state_[j];
    if (c.reserved) {
        throw std::runtime_error(
            "truesum: a device column may only be reserved once");
    }

    c.exponent = exponent;
    c.nlimbs = limbs_for_bits(bits) < 1 ? 1 : limbs_for_bits(bits);
    c.reserved = true;

    if (0 != c.rows) {
        // One allocation per limb position. The stream-ordered pool is what
        // makes that affordable -- plain cudaMalloc is tens of microseconds a
        // call, and a wide column needs tens of them. Widening later appends
        // to `bases` without disturbing anything already allocated.
        const std::size_t bytes = c.rows * sizeof(limb_t);
        c.bases.resize(c.nlimbs);
        for (std::size_t k = 0; k < c.nlimbs; ++k) {
            // On this accumulation matrix's stream, not the null stream:
            // stream-ordered memory is valid for work ordered after the
            // allocation *in that stream*, and with a caller-supplied
            // non-blocking stream nothing orders the null stream against it.
            cuda(MallocAsync(&c.bases[k], bytes, st_compute_));
            // Zeroed later, all of them together -- see flush_pending_zero.
            zero_ptr_.push_back(c.bases[k]);
            zero_rows_.push_back(c.rows);
        }
        cuda(Malloc(&c.dev_bases, c.nlimbs * sizeof(limb_t *)));
        cuda(Memcpy(c.dev_bases, c.bases.data(), c.nlimbs * sizeof(limb_t *),
                    cudaMemcpyHostToDevice));
    }
    descriptors_stale_ = true;
}

// Appends limb positions until the column has `needed` of them. Existing
// arrays are not touched -- that is what limb-major buys, and what the
// per-limb-position allocation preserves -- so growth is an allocation plus a
// sign fill, with no data movement at all.
//
// Synchronizes first. Widening rewrites the descriptor array that in-flight
// kernels are reading, and it is rare enough that draining once is simpler
// than versioning the descriptors.
void
CudaAccumulationMatrix::grow_column(std::size_t j, std::size_t needed)
{
    Column &c = cols_state_[j];
    if (c.nlimbs >= needed) return;
    flush_pending_zero();
    synchronize();

    const std::size_t bytes = c.rows * sizeof(limb_t);
    const dim3 grid = launch_grid(c.rows, 1);
    for (std::size_t k = c.nlimbs; k < needed; ++k) {
        limb_t *fresh = nullptr;
        cuda(MallocAsync(&fresh, bytes, st_compute_));
        sign_fill_kernel<<<grid, kBlock, 0, st_compute_>>>(
            fresh, c.bases[k - 1], c.rows);
        c.bases.push_back(fresh);
    }
    c.nlimbs = needed;

    // A fresh array rather than an overwrite: the old one may still be under a
    // kernel that has not retired, and cudaFreeAsync releases it in order.
    limb_t **fresh_bases = nullptr;
    cuda(Malloc(&fresh_bases, c.nlimbs * sizeof(limb_t *)));
    cuda(Memcpy(fresh_bases, c.bases.data(), c.nlimbs * sizeof(limb_t *),
                cudaMemcpyHostToDevice));
    cudaFree(c.dev_bases);
    c.dev_bases = fresh_bases;
    descriptors_stale_ = true;
}

// Lowers a column's exponent, shifting every entry left to match. Mirrors
// AccumulationMatrix::rescale: widen enough to hold the shifted value, then
// shift. Unlike widening this moves every bit in the column, which is why the
// exponent is the half worth pre-sizing correctly.
void
CudaAccumulationMatrix::rescale_column(std::size_t j, int new_exponent)
{
    Column &c = cols_state_[j];
    if (new_exponent >= c.exponent) return;
    const unsigned shift = static_cast<unsigned>(
        static_cast<long long>(c.exponent) - new_exponent);

    const std::size_t widened =
        c.max_addend_bits + shift + ceil_log2(c.add_count + 1) + 1;
    const std::size_t ndst =
        std::max(limbs_for_bits(widened), c.nlimbs + shift / limbs::kLimbBits);
    grow_column(j, ndst);  // synchronizes, and sign-fills what it appends

    shift_left_kernel<<<launch_grid(c.rows, 1), kBlock, 0, st_compute_>>>(
        c.dev_bases, static_cast<unsigned>(c.nlimbs), c.rows, shift);

    c.max_addend_bits += shift;
    c.exponent = new_exponent;
    // The exponent lives in the device descriptor, and grow_column only marks
    // it stale when it actually appends. A rescale that needed no new limbs
    // would otherwise leave the kernel computing shifts against the old
    // exponent -- which goes negative, makes `off` enormous, and drops every
    // value in the batch without a word.
    descriptors_stale_ = true;
}

// The derived bound AccumulationMatrix::fit_column applies, checked against
// what the column was actually reserved for. Called once per column per batch,
// after the survey has established that batch's extent.
void
CudaAccumulationMatrix::require_fit(std::size_t j, int min_exponent,
                                    int max_top, std::size_t count)
{
    if (min_exponent < cols_state_[j].exponent) {
        rescale_column(j, static_cast<int>(min_exponent));
    }
    Column &c = cols_state_[j];
    c.max_addend_bits = std::max(
        c.max_addend_bits, static_cast<std::size_t>(max_top - c.exponent));
    c.add_count += count;

    const std::size_t needed_bits =
        c.max_addend_bits + ceil_log2(c.add_count + 1) + 1;
    const std::size_t needed = limbs_for_bits(needed_bits);
    if (needed > c.nlimbs) grow_column(j, needed);
}

// Shared tail of both reserve_for entry points: turn per-column exponent
// extents into reservations. Mirrors AccumulationMatrix::reserve_for, with one
// deliberate difference -- a column with nothing in it is still reserved,
// since the device cannot grow one later.
void
CudaAccumulationMatrix::reserve_from_extents(const int *low, const int *high,
                                             const char *any, std::size_t count)
{
    const std::size_t headroom =
        ceil_log2(count < 1 ? 1 : count) + 1;  // +1 for the sign
    for (std::size_t j = 0; j < cols_; ++j) {
        if (!any[j]) {
            reserve_column(j, 0, limbs::kLimbBits);
            continue;
        }
        reserve_column(j, low[j],
                       static_cast<std::size_t>(high[j] - low[j]) + headroom);
    }
}

void
CudaAccumulationMatrix::reserve_for(const double *b, std::size_t count,
                                    std::size_t col_stride)
{
    if (0 == cols_) return;
    const std::size_t stride = col_stride ? col_stride : rows_;
    std::vector<int> low(cols_, 0), high(cols_, 0);
    std::vector<char> any(cols_, 0);

    for (std::size_t j = 0; j < cols_; ++j) {
        // A column with no scale yet needs the exact true-ulp minimum, so the
        // cutoff is one nothing can clear.
        const Column &cj = cols_state_[j];
        const kernels::Survey sv =
            kernels::survey()(b + j * stride + cj.first_row, cj.rows,
                              std::numeric_limits<long long>::max());
        if (sv.nonfinite) {
            throw std::domain_error(
                "truesum: cannot reserve from a non-finite value");
        }
        if (!sv.any) continue;
        low[j] = sv.min_exponent;
        high[j] = sv.max_top;
        any[j] = 1;
    }
    reserve_from_extents(low.data(), high.data(), any.data(), count);
}

// Writing the sentinels with a kernel rather than copying them up. The copy
// this replaces was pageable and host-to-device, which the driver is allowed to
// stage through a pinned buffer of its own, and staging one may synchronize the
// stream it was queued on -- a drain in the middle of an accumulate, to deliver
// 768 bytes at 64 columns.
void
CudaAccumulationMatrix::begin_survey()
{
    const unsigned n = static_cast<unsigned>(cols_);
    const unsigned blocks = (n + kBlock - 1) / kBlock;
    survey_init_kernel<<<blocks, kBlock, 0, st_compute_>>>(
        static_cast<Survey *>(survey_out_), n);
}

// The verdict comes back into pinned memory, so this copy is a real DMA. The
// synchronize is not: the host cannot size a column until it has read the
// extents, and that is the one wait on this path that nothing can remove.
const Survey *
CudaAccumulationMatrix::end_survey()
{
    cuda(MemcpyAsync(survey_host_, survey_out_, cols_ * sizeof(Survey),
                     cudaMemcpyDeviceToHost, st_compute_));
    synchronize();
    return static_cast<const Survey *>(survey_host_);
}

void
CudaAccumulationMatrix::reserve_for_device(const double *b, std::size_t count,
                                           std::size_t col_stride)
{
    if (0 == cols_ || 0 == rows_) return;
    const std::size_t stride = col_stride ? col_stride : rows_;

    begin_survey();
    survey_kernel<<<launch_grid(rows_, cols_), kBlock, 0, st_compute_>>>(
        b, static_cast<const ColumnShape *>(shapes_), 0, stride,
        static_cast<Survey *>(survey_out_));
    const Survey *surveys = end_survey();

    std::vector<int> low(cols_, 0), high(cols_, 0);
    std::vector<char> any(cols_, 0);
    for (std::size_t j = 0; j < cols_; ++j) {
        if (surveys[j].nonfinite) {
            throw std::domain_error(
                "truesum: cannot reserve from a non-finite value");
        }
        if (!surveys[j].any) continue;
        low[j] = surveys[j].min_exponent;
        high[j] = surveys[j].max_top;
        any[j] = 1;
    }
    reserve_from_extents(low.data(), high.data(), any.data(), count);
}

void
CudaAccumulationMatrix::reserve_like(const AccumulationMatrix &cpu)
{
    if (cpu.rows() != rows_ || cpu.cols() != cols_) {
        throw std::runtime_error(
            "truesum: reserve_like requires matching dimensions");
    }
    // Same shape, not merely the same extent: a triangular column is a
    // different length, so copying a full accumulation matrix's widths into a
    // triangular one would reserve the right bits for the wrong entries.
    if (cpu.symmetric() != symmetric_ ||
        (symmetric_ && cpu.uplo() != uplo_)) {
        throw std::runtime_error(
            "truesum: reserve_like requires the same symmetry and uplo");
    }
    for (std::size_t j = 0; j < cols_; ++j) {
        reserve_column(j, cpu.column_exponent(j), cpu.column_bit_width(j));
    }
}

int
CudaAccumulationMatrix::column_occupancy(std::size_t j) const
{
    if (j >= cols_)
        throw std::out_of_range("truesum: column index out of range");
    flush_pending_zero();
    synchronize();
    const_cast<CudaAccumulationMatrix *>(this)->harvest_occupancy();
    return cols_state_[j].max_limb_used;
}

// The device staging accumulates across launches, so the mapped array already
// holds the high-water mark; this only copies it where the column keeps it.
void
CudaAccumulationMatrix::harvest_occupancy()
{
    if (nullptr == occupancy_host_) return;
    for (std::size_t j = 0; j < cols_; ++j) {
        cols_state_[j].max_limb_used = occupancy_host_[j];
    }
}

std::size_t
CudaAccumulationMatrix::memory_bytes() const
{
    std::size_t total = 0;
    for (const auto &c : cols_state_) {
        total += c.nlimbs * c.rows * sizeof(limb_t);
    }
    return total;
}

// One launch over every array still owing a zero. The y extent is one block
// row per array, so a short column costs a block that exits at once rather
// than a separate API call that does not.
void
CudaAccumulationMatrix::flush_pending_zero() const
{
    if (zero_ptr_.empty()) return;
    const std::size_t n = zero_ptr_.size();

    if (zero_capacity_ < n) {
        cudaFree(zero_targets_);
        zero_targets_ = nullptr;
        cuda(Malloc(&zero_targets_, n * sizeof(ZeroTarget)));
        zero_capacity_ = n;
    }
    std::vector<ZeroTarget> host(n);
    std::size_t widest = 0;
    for (std::size_t i = 0; i < n; ++i) {
        host[i].p = zero_ptr_[i];
        host[i].rows = static_cast<unsigned>(zero_rows_[i]);
        if (zero_rows_[i] > widest) widest = zero_rows_[i];
    }
    // Synchronous, and on the null stream rather than st_compute_, which is
    // deliberate on all three counts. The host array is pageable, so an async
    // copy could stall the stream to stage it; being synchronous, this one is
    // complete before the launch below reads it, which is the ordering that
    // matters; and the cudaFree/cudaMalloc a few lines above synchronize far
    // harder anyway. Nothing here runs unless zeroing is actually pending.
    cuda(Memcpy(zero_targets_, host.data(), n * sizeof(ZeroTarget),
                cudaMemcpyHostToDevice));

    unsigned gx = static_cast<unsigned>((widest + kBlock - 1) / kBlock);
    if (gx > kMaxBlocks) gx = kMaxBlocks;
    if (0 == gx) gx = 1;
    zero_kernel<<<dim3(gx, static_cast<unsigned>(n)), kBlock, 0, st_compute_>>>(
        static_cast<const ZeroTarget *>(zero_targets_));

    zero_ptr_.clear();
    zero_rows_.clear();
}

void
CudaAccumulationMatrix::sync_descriptors()
{
    if (!descriptors_stale_) return;

    // The staging buffer may still be under a previous upload. In practice it
    // never is -- descriptors only go stale from grow_column and
    // rescale_column, both of which have already drained the stream -- but the
    // guard costs nothing when the event has long since fired.
    if (desc_in_flight_) {
        cuda(EventSynchronize(ev_desc_));
        desc_in_flight_ = false;
    }
    ColumnDesc *host = static_cast<ColumnDesc *>(desc_host_);
    for (std::size_t j = 0; j < cols_; ++j) {
        host[j].bases = cols_state_[j].dev_bases;
        host[j].exponent = cols_state_[j].exponent;
        host[j].nlimbs = static_cast<unsigned>(cols_state_[j].nlimbs);
    }
    // On st_compute_, so it is ordered before the launch that reads it without
    // the host waiting for anything already queued there.
    cuda(MemcpyAsync(descriptors_, host, cols_ * sizeof(ColumnDesc),
                     cudaMemcpyHostToDevice, st_compute_));
    cuda(EventRecord(ev_desc_, st_compute_));
    desc_in_flight_ = true;
    descriptors_stale_ = false;
}

void
CudaAccumulationMatrix::ensure_slots(std::size_t words)
{
    if (slot_words_ >= words) return;
    synchronize();
    for (Slot &s : slots_) {
        if (s.mapped) cudaFreeHost(s.mapped);
        s.mapped = nullptr;
        cuda(HostAlloc(&s.mapped, words * sizeof(double), cudaHostAllocMapped));
        if (!s.ev_done) {
            cuda(EventCreateWithFlags(&s.ev_done, cudaEventDisableTiming));
        }
        s.in_flight = false;
    }
    slot_words_ = words;
}

// Validates a packed column-major batch on the host, using the same survey the
// CPU accumulation matrix uses -- which is the AVX-512 one where available, at
// ~0.12 ns/elem. Cheap enough to hide entirely behind the transfer it runs
// against.
void
CudaAccumulationMatrix::require_all_reserved() const
{
    for (std::size_t j = 0; j < cols_; ++j) {
        if (!cols_state_[j].reserved) {
            throw std::runtime_error(
                "truesum: every device column must be reserved before "
                "accumulation");
        }
    }
}

void
CudaAccumulationMatrix::validate_host_survey(const double *b,
                                             std::size_t col_stride,
                                             const Survey *supplied)
{
    require_all_reserved();

    // Not threaded, and measured rather than assumed: with two slots in
    // rotation the survey of the next batch already overlaps the kernel
    // streaming the previous one -- 0.5 ms of survey inside 1.26 ms of
    // transfer at 65536x64 -- so it is not on the critical path. Threading it
    // measured flat on this path, helped one shape of the staged path and hurt
    // another.
    for (std::size_t j = 0; j < cols_; ++j) {
        const Column &cj = cols_state_[j];
        const kernels::Survey sv =
            nullptr != supplied
                ? supplied[j]
                : kernels::survey()(b + j * col_stride + cj.first_row, cj.rows,
                                    cj.exponent);
        if (sv.nonfinite) {
            throw std::domain_error(
                "truesum: cannot accumulate a non-finite value");
        }
        if (!sv.any) continue;
        require_fit(j, sv.min_exponent, sv.max_top);
    }
}

// The host path runs the transfer and the survey against each other: stage into
// pinned memory, start the copy, then survey that same buffer while the DMA is
// in flight. The survey is pure host work on host memory, so it costs nothing
// the transfer was not already going to spend.
//
// Validation therefore still happens *before* the accumulate is launched, so a
// bad batch is rejected without having touched the accumulation matrix -- which
// a device-side survey cannot do without a round-trip that drains the pipeline.

double *
CudaAccumulationMatrix::acquire_input()
{
    ensure_slots(rows_ * cols_);
    Slot &s = slots_[slot_];
    if (s.in_flight) {
        cuda(EventSynchronize(s.ev_done));
        s.in_flight = false;
    }
    return s.mapped;
}

void
CudaAccumulationMatrix::add_matrix_col_major(const double *b,
                                             const Survey *surveys,
                                             std::size_t col_stride)
{
    if (0 == rows_ || 0 == cols_) return;
    const std::size_t stride = col_stride ? col_stride : rows_;

    // The kernel dereferences this on the device, so pageable memory faults
    // rather than merely running slowly. Checked here, where the diagnostic
    // can say what is wrong and what to allocate instead -- and checked at
    // all because this path no longer has a staging copy to fall back on.
    cudaPointerAttributes attr{};
    const cudaError_t st = cudaPointerGetAttributes(&attr, b);
    if (cudaSuccess != st || nullptr == attr.devicePointer) {
        cudaGetLastError();  // an unregistered pointer leaves this sticky
        throw std::invalid_argument(
            "truesum: add_matrix_col_major needs page-locked, "
            "device-mapped host memory -- cudaHostAlloc with "
            "cudaHostAllocMapped, or cudaHostRegister with "
            "cudaHostRegisterMapped");
    }

    if (presized_) {
        if (++submitted_matrices_ > declared_matrices_) {
            throw std::runtime_error(
                "truesum: more matrices accumulated than were described to the "
                "constructor; the width bound holds for that many and no more");
        }
    } else {
        // Surveyed through the host pointer; launched through the device one,
        // which is the same address under unified addressing but need not be
        // for a registered range.
        validate_host_survey(b, stride, surveys);
    }
    { const double *one = static_cast<const double *>(attr.devicePointer);
      launch_accumulate(&one, 1, stride); }

    // A buffer from acquire_input is one of ours, and its slot's event is what
    // that call blocks on. Keep the rotation's bookkeeping current so both
    // ways of waiting stay correct: the handle below, and acquire_input.
    Slot &s = slots_[slot_];
    if (b == s.mapped && rows_ == stride) {
        cuda(EventRecord(s.ev_done, st_compute_));
        s.in_flight = true;
        slot_ ^= 1;
    }

}


// One survey staging, one launch per matrix. The kernel reduces with atomicMin
// and atomicMax, so successive launches accumulate the union of their extents
// -- which is exactly the aggregate a batch has to be sized for.
void
CudaAccumulationMatrix::survey_device_inputs(const double *const *b,
                                             std::size_t count,
                                             std::size_t col_stride,
                                             const Survey *const *supplied)
{
    if (nullptr != supplied) {
        // Aggregated on the host from what the caller already knows, so no
        // kernel is launched and nothing is waited for.
        for (std::size_t j = 0; j < cols_; ++j) {
            for (std::size_t k = 0; k < count; ++k) {
                if (supplied[k][j].nonfinite) {
                    throw std::domain_error(
                        "truesum: cannot accumulate a non-finite value");
                }
                if (!supplied[k][j].any) continue;
                require_fit(j, supplied[k][j].min_exponent,
                            supplied[k][j].max_top, 1);
            }
        }
        return;
    }
    const dim3 grid = launch_grid(rows_, cols_);
    begin_survey();
    for (std::size_t k = 0; k < count; ++k) {
        survey_kernel<<<grid, kBlock, 0, st_compute_>>>(
            b[k], static_cast<const ColumnShape *>(shapes_), 0, col_stride,
            static_cast<Survey *>(survey_out_));
    }
    const Survey *surveys = end_survey();

    for (std::size_t j = 0; j < cols_; ++j) {
        if (surveys[j].nonfinite) {
            throw std::domain_error(
                "truesum: cannot accumulate a non-finite value");
        }
        if (!surveys[j].any) continue;
        require_fit(j, surveys[j].min_exponent, surveys[j].max_top, count);
    }
}

void
CudaAccumulationMatrix::add_matrices_col_major_device(const double *const *b,
                                                      std::size_t count,
                                                      const Survey *const *surveys,
                                                      std::size_t col_stride)
{
    if (0 == count || 0 == rows_ || 0 == cols_) return;
    if (count > kMaxBatchInputs) {
        std::ostringstream os;
        os << "truesum: at most " << kMaxBatchInputs
           << " matrices may be batched into one pass; " << count << " given";
        throw std::invalid_argument(os.str());
    }
    require_all_reserved();
    const std::size_t stride = col_stride ? col_stride : rows_;

    if (presized_) {
        submitted_matrices_ += count;
        if (submitted_matrices_ > declared_matrices_) {
            throw std::runtime_error(
                "truesum: more matrices accumulated than were described to the "
                "constructor; the width bound holds for that many and no more");
        }
    } else {
        // Every matrix must be surveyed and the column fitted before any of
        // them is added: a widen partway through the batch would leave earlier
        // matrices already written into a column of the wrong shape.
        survey_device_inputs(b, count, stride, surveys);
    }
    sync_descriptors();
    launch_accumulate(b, count, stride);
}

void
CudaAccumulationMatrix::add_matrix_col_major_device(const double *b,
                                                    const Survey *surveys,
                                                    std::size_t col_stride)
{
    if (0 == rows_ || 0 == cols_) return;
    accumulate_device(b, col_stride ? col_stride : rows_, surveys);
}

CUstream_st *
CudaAccumulationMatrix::stream() const
{
    return st_compute_;
}

void
CudaAccumulationMatrix::synchronize() const
{
    cuda(StreamSynchronize(st_compute_));
    report_contradictions();
}

unsigned
CudaAccumulationMatrix::column_contradictions(std::size_t j) const
{
    if (j >= cols_)
        throw std::out_of_range("truesum: column index out of range");
    cuda(StreamSynchronize(st_compute_));  // not synchronize(): do not throw
    return nullptr == flags_host_ ? 0u : flags_host_[j];
}

// A batch reached outside what its column was sized for. Pre-sized from
// producer metadata that means the metadata was wrong; sized from a survey
// taken here it means the survey and the accumulate disagree, which is a
// bug. Either way the sums are already wrong, so this reports rather than
// recovers -- the same contract, and the same wording, as the CPU container.
void
CudaAccumulationMatrix::report_contradictions() const
{
    if (nullptr == flags_host_) return;
    for (std::size_t j = 0; j < cols_; ++j) {
        const unsigned flags = flags_host_[j];
        if (0 == flags) continue;
        std::ostringstream os;
        os << "truesum: column " << j
           << " received values it was not sized for:";
        if (0 != (flags & kBadNonFinite)) os << " a non-finite value;";
        if (0 != (flags & kBadExponent)) {
            os << " an exponent below the column's;";
        }
        if (0 != (flags & kBadWidth)) os << " an addend past its width;";
        if (0 != (flags & kBadOverflow)) os << " a sum past its width;";
        os << " the accumulation matrix is no longer consistent";
        throw std::runtime_error(os.str());
    }
}

void
CudaAccumulationMatrix::accumulate_device(const double *b,
                                          std::size_t col_stride,
                                          const Survey *supplied)
{
    for (std::size_t j = 0; j < cols_; ++j) {
        if (!cols_state_[j].reserved) {
            throw std::runtime_error(
                "truesum: every device column must be reserved before "
                "accumulation");
        }
    }
    if (presized_) {
        // The drain below is this path's fixed cost, and pre-sizing is what
        // removes it: with the extents known there is no question to ask the
        // device, so the accumulate launches straight away and this call
        // never blocks.
        if (++submitted_matrices_ > declared_matrices_) {
            throw std::runtime_error(
                "truesum: more matrices accumulated than were described to "
                "the constructor; the width bound holds for that many and "
                "no more");
        }
        launch_accumulate(&b, 1, col_stride);
        return;
    }
    sync_descriptors();

    const dim3 grid = launch_grid(rows_, cols_);

    // First pass: learn each column's exponent range, and reject anything the
    // reservation cannot hold. The decisions the CPU makes by rescaling and
    // widening are errors here, because neither is possible mid-launch.
    // A supplied survey removes both the kernel and the wait: the extents are
    // already on the host, so there is nothing to ask the device and nothing to
    // come back. That drain is this path's fixed cost, 19.9 us a batch.
    const Survey *surveys = supplied;
    if (nullptr == surveys) {
        begin_survey();
        survey_kernel<<<grid, kBlock, 0, st_compute_>>>(
            b, static_cast<const ColumnShape *>(shapes_), 0, col_stride,
            static_cast<Survey *>(survey_out_));
        // With the input already on the device, there is nothing to survey on
        // the host, so the verdict has to come back before the accumulate can
        // be allowed to run.
        surveys = end_survey();
    }

    for (std::size_t j = 0; j < cols_; ++j) {
        if (surveys[j].nonfinite) {
            throw std::domain_error(
                "truesum: cannot accumulate a non-finite value");
        }
        if (!surveys[j].any) continue;
        require_fit(j, surveys[j].min_exponent, surveys[j].max_top);
    }

    launch_accumulate(&b, 1, col_stride);
}

// Queues the accumulate. Asynchronous: the caller is not blocked, so batches
// pipeline against each other without the caller doing anything.
void
CudaAccumulationMatrix::launch_accumulate(const double *const *b,
                                          std::size_t count,
                                          std::size_t col_stride)
{
    flush_pending_zero();
    for (std::size_t j = 0; j < cols_; ++j) {
        if (!cols_state_[j].reserved) {
            throw std::runtime_error(
                "truesum: every device column must be reserved before "
                "accumulation");
        }
    }
    sync_descriptors();

    const dim3 grid = launch_grid(rows_, cols_);

    InputSet in;
    in.n = static_cast<unsigned>(count);
    for (std::size_t k = 0; k < count; ++k) in.p[k] = b[k];
    accumulate_kernel<64><<<grid, kBlock, 0, st_compute_>>>(
        static_cast<const ColumnDesc *>(descriptors_), in,
        static_cast<const ColumnShape *>(shapes_), col_stride,
        occupancy_device_, occupancy_host_, flags_device_, flags_host_,
        ticket_, static_cast<unsigned>(cols_));
}

std::vector<limb_t>
CudaAccumulationMatrix::entry_limbs(std::size_t i, std::size_t j) const
{
    check_index(i, j);
    flush_pending_zero();
    synchronize();
    std::size_t col = 0, slot = 0;
    locate(col, slot, i, j);
    const Column &c = cols_state_[col];
    std::vector<limb_t> v(c.nlimbs);
    // One small copy per limb position. This is a readback path, not a hot
    // one; a bulk download would gather whole limb arrays instead.
    for (std::size_t k = 0; k < c.nlimbs; ++k) {
        cuda(Memcpy(&v[k], c.bases[k] + slot, sizeof(limb_t),
                    cudaMemcpyDeviceToHost));
    }
    return v;
}

std::vector<limb_t>
CudaAccumulationMatrix::download_column(std::size_t j) const
{
    if (j >= cols_)
        throw std::out_of_range("truesum: column index out of range");
    flush_pending_zero();
    synchronize();
    const Column &c = cols_state_[j];
    std::vector<limb_t> out(c.nlimbs * c.rows);
    for (std::size_t k = 0; k < c.nlimbs; ++k) {
        cuda(Memcpy(out.data() + k * c.rows, c.bases[k],
                    c.rows * sizeof(limb_t), cudaMemcpyDeviceToHost));
    }
    return out;
}

}  // namespace truesum

#undef cuda
