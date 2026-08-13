#include <cuda_runtime.h>

#include <climits>
#include <cstring>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "cbfp/column_accumulator.hpp"
#include "cbfp/cuda_accumulator.hpp"
#include "kernels.hpp"

namespace cbfp {
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
    os << "cbfp: " << call << " failed at " << file << ":" << line << ": "
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
// diagnostic, so a failure names the function rather than an opaque
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

// Blocks in the whole grid, not just its x extent. Both kernels grid-stride,
// so any smaller grid stays correct -- it only gives each thread more rows.
//
// Capping the total is what makes the survey's reduction pay for itself. Left
// uncapped, 65536 rows over 64 columns launches 16384 blocks, each folding a
// single value per thread and then paying a full eight-step shared-memory tree
// to do it. Capped, each thread folds sixteen rows serially first and the tree
// is amortised across them: measured 195 -> 113 us at that shape, and
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

// What the first pass learns about a column, mirroring kernels::Survey. Reduced
// across the whole column with atomics, so every field is an atomic-friendly
// type rather than the host struct's bools.
//
// The exponents are int where the host struct uses long long. A double's
// exponent lives in [-1074, 1077] and the container caps its own at 2^24, so
// 32 bits is three orders of magnitude more than enough -- and it buys 32-bit
// atomics and halves this kernel's shared memory.
struct DeviceSurvey {
    int min_exponent;
    int max_top;
    int any;
    int nonfinite;
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

// The device twin of the scalar kernel's `split`, field for field. Kept in
// lockstep with it deliberately: the two implementations have to agree about
// what a double means before they can be compared bit for bit.
//
//   normal    (biased != 0): m = 2^52 | frac, e = biased - 1075
//   subnormal (biased == 0): m = frac,        e = -1074
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
// a column's scale is settled. That does not pay here: this pass is bound by
// reading the column, and the trailing-zero count is a few ALU ops on a value
// already in registers. One exact pass is both simpler and cheaper.
__global__ void
survey_kernel(const double *__restrict__ values, std::size_t rows,
              std::size_t col_stride, DeviceSurvey *__restrict__ out)
{
    __shared__ int s_min[kBlock];
    __shared__ int s_max[kBlock];
    __shared__ int s_any[kBlock];
    __shared__ int s_bad[kBlock];

    const unsigned tid = threadIdx.x;
    const unsigned j = blockIdx.y;
    const double *col = values + static_cast<std::size_t>(j) * col_stride;

    int tmin = INT_MAX;
    int tmax = INT_MIN;
    int tany = 0;
    int tbad = 0;

    const std::size_t step = static_cast<std::size_t>(gridDim.x) * blockDim.x;
    for (std::size_t i =
             static_cast<std::size_t>(blockIdx.x) * blockDim.x + tid;
         i < rows; i += step) {
        unsigned long long m;
        int e, top;
        bool neg, bad;
        split_device(col[i], m, e, top, neg, bad);
        if (bad) tbad = 1;
        if (0 != m) {
            tany = 1;
            if (e < tmin) tmin = e;
            if (top > tmax) tmax = top;
        }
    }

    s_min[tid] = tmin;
    s_max[tid] = tmax;
    s_any[tid] = tany;
    s_bad[tid] = tbad;
    __syncthreads();
    for (unsigned s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            if (s_min[tid + s] < s_min[tid]) s_min[tid] = s_min[tid + s];
            if (s_max[tid + s] > s_max[tid]) s_max[tid] = s_max[tid + s];
            s_any[tid] |= s_any[tid + s];
            s_bad[tid] |= s_bad[tid + s];
        }
        __syncthreads();
    }

    // One atomic per block rather than one per thread.
    if (0 == tid) {
        if (s_any[0]) {
            atomicMin(&out[j].min_exponent, s_min[0]);
            atomicMax(&out[j].max_top, s_max[0]);
            atomicOr(&out[j].any, 1);
        }
        if (s_bad[0]) atomicOr(&out[j].nonfinite, 1);
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
accumulate_kernel(const ColumnDesc *__restrict__ cols,
                  const double *__restrict__ values, std::size_t rows,
                  std::size_t col_stride, int *__restrict__ occupancy_device,
                  int *__restrict__ occupancy_host,
                  unsigned *__restrict__ ticket, unsigned ncols)
{
    // The addend split, the offset arithmetic and the limb mask below are all
    // written in terms of kRadix. What is not yet written is the carry-save
    // apply: at a reduced radix the inner loop becomes a bare add with no
    // carry-out test and no dependency between limb positions, plus a
    // normalization pass before readback. That is the next milestone.
    static_assert(kRadix == 64,
                  "only the canonical radix is implemented; carry-save at "
                  "radix 52 is the next milestone");
    constexpr unsigned long long kLimbMask =
        kRadix == 64 ? ~0ull : ((1ull << kRadix) - 1);

    const unsigned j = blockIdx.y;
    const ColumnDesc c = cols[j];
    const double *col = values + static_cast<std::size_t>(j) * col_stride;

    __shared__ int s_lim[kBlock];
    const unsigned tid = threadIdx.x;
    // Highest limb position this thread disturbs. The carry loop already knows
    // where it stopped, so this costs a comparison rather than a pass.
    int t_lim = -1;

    const std::size_t step = static_cast<std::size_t>(gridDim.x) * blockDim.x;
    for (std::size_t i =
             static_cast<std::size_t>(blockIdx.x) * blockDim.x + tid;
         i < rows; i += step) {
        unsigned long long m;
        int e, top;
        bool neg, bad;
        split_device(col[i], m, e, top, neg, bad);
        if (0 == m) continue;

        // The host has already checked that every value fits the column it
        // was reserved for, so the shift cannot be negative. Were it ever
        // negative anyway, the unsigned conversion puts `off` far above
        // nlimbs and the loop below simply does not run.
        const int shift = e - c.exponent;
        const unsigned off = static_cast<unsigned>(shift) / kRadix;
        const unsigned bit = static_cast<unsigned>(shift) % kRadix;

        // A 53-bit significand at intra-limb offset `bit` spans two limbs at
        // radix 64 and at radix 52 alike.
        unsigned long long lo, hi;
        if (0 == bit) {
            lo = m & kLimbMask;
            hi = kRadix == 64 ? 0ull : (m >> kRadix) & kLimbMask;
        } else {
            lo = (m << bit) & kLimbMask;
            hi = (m >> (kRadix - bit)) & kLimbMask;
        }

        // Add and subtract share the loop shape; `carry` is a borrow when the
        // lane is negative. Past the addend with nothing propagating, no
        // higher limb can change.
        unsigned long long carry = 0;
        for (unsigned p = off; p < c.nlimbs; ++p) {
            const unsigned long long a =
                (p == off) ? lo : ((p == off + 1) ? hi : 0ull);
            limb_t *dst = c.bases[p] + i;
            const unsigned long long x = *dst;
            unsigned long long written;
            if (neg) {
                const unsigned long long d = x - a;
                const unsigned long long b1 = (x < a) ? 1ull : 0ull;
                const unsigned long long d2 = d - carry;
                const unsigned long long b2 = (d < carry) ? 1ull : 0ull;
                written = d2;
                carry = b1 | b2;
            } else {
                const unsigned long long s = x + a;
                const unsigned long long c1 = (s < x) ? 1ull : 0ull;
                const unsigned long long s2 = s + carry;
                const unsigned long long c2 = (s2 < s) ? 1ull : 0ull;
                written = s2;
                carry = c1 | c2;
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

    // Fold the per-thread maxima, reduce across blocks in device memory with
    // ordinary atomics, then elect one block to carry the finished array to
    // host memory. Atomics never cross PCIe: measured, that costs 260x.
    s_lim[tid] = t_lim;
    __syncthreads();
    for (unsigned s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s && s_lim[tid + s] > s_lim[tid]) s_lim[tid] = s_lim[tid + s];
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
    }
    if (0 == tid) *ticket = 0;
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

}  // namespace

bool
cuda_available()
{
    int n = 0;
    return cudaSuccess == cudaGetDeviceCount(&n) && n > 0;
}

CudaColumnBlockMatrix::CudaColumnBlockMatrix(std::size_t rows, std::size_t cols)
    : rows_(rows), cols_(cols), cols_state_(cols)
{
    if (!cuda_available()) {
        throw std::runtime_error("cbfp: no usable CUDA device");
    }
    // Mapping has to be enabled before the context exists, so this fails
    // harmlessly if one is already active with the flag set. The cudaHostAlloc
    // below is the check that matters -- it fails loudly if mapping is really
    // unavailable.
    cudaSetDeviceFlags(cudaDeviceMapHost);
    // Columns become the grid's y dimension, the only launch parameter here
    // that is not fixed at compile time or clamped. Bounding it once, by name,
    // is what lets every launch below go unchecked.
    int max_grid_y = 0;
    cuda(DeviceGetAttribute(&max_grid_y, cudaDevAttrMaxGridDimY, 0));
    if (cols_ > static_cast<std::size_t>(max_grid_y)) {
        std::ostringstream os;
        os << "cbfp: " << cols_ << " columns exceeds this device's grid y "
           << "limit of " << max_grid_y;
        throw std::runtime_error(os.str());
    }
    cuda(StreamCreate(&st_copy_));
    cuda(StreamCreate(&st_compute_));

    if (0 != cols_) {
        cuda(Malloc(&occupancy_device_, cols_ * sizeof(int)));
        cuda(Memset(occupancy_device_, 0xFF, cols_ * sizeof(int)));  // -1
        cuda(HostAlloc(&occupancy_host_, cols_ * sizeof(int),
                       cudaHostAllocMapped));
        for (std::size_t j = 0; j < cols_; ++j) occupancy_host_[j] = -1;
        cuda(Malloc(&ticket_, sizeof(unsigned)));
        cuda(Memset(ticket_, 0, sizeof(unsigned)));
        cuda(Malloc(&descriptors_, cols_ * sizeof(ColumnDesc)));
        cuda(Malloc(&survey_out_, cols_ * sizeof(DeviceSurvey)));
    }
}

CudaColumnBlockMatrix::~CudaColumnBlockMatrix()
{
    // Deliberately unchecked: a destructor must not throw, and there is
    // nothing useful to do about a failed free during teardown.
    for (auto &c : cols_state_) {
        for (limb_t *p : c.bases) cudaFreeAsync(p, 0);
        cudaFree(c.dev_bases);
    }
    cudaFree(descriptors_);
    cudaFree(survey_out_);
    cudaFree(occupancy_device_);
    cudaFreeHost(occupancy_host_);
    cudaFree(ticket_);
    for (Slot &s : slots_) {
        if (s.pinned) cudaFreeHost(s.pinned);
        if (s.device) cudaFree(s.device);
        if (s.ev_copied) cudaEventDestroy(s.ev_copied);
        if (s.ev_done) cudaEventDestroy(s.ev_done);
    }
    if (st_copy_) cudaStreamDestroy(st_copy_);
    if (st_compute_) {
        cudaStreamDestroy(st_compute_);
    }
}

void
CudaColumnBlockMatrix::check_index(std::size_t i, std::size_t j) const
{
    if (i >= rows_ || j >= cols_) {
        throw std::out_of_range("cbfp: matrix index out of range");
    }
}

int
CudaColumnBlockMatrix::column_exponent(std::size_t j) const
{
    if (j >= cols_) throw std::out_of_range("cbfp: column index out of range");
    return cols_state_[j].exponent;
}

std::size_t
CudaColumnBlockMatrix::column_limbs(std::size_t j) const
{
    if (j >= cols_) throw std::out_of_range("cbfp: column index out of range");
    return cols_state_[j].nlimbs;
}

void
CudaColumnBlockMatrix::reserve_column(std::size_t j, int exponent,
                                      std::size_t bits)
{
    if (j >= cols_) throw std::out_of_range("cbfp: column index out of range");
    Column &c = cols_state_[j];
    if (c.reserved) {
        throw std::runtime_error(
            "cbfp: a device column may only be reserved once");
    }

    c.exponent = exponent;
    c.nlimbs = limbs_for_bits(bits) < 1 ? 1 : limbs_for_bits(bits);
    c.reserved = true;

    if (0 != rows_) {
        // One allocation per limb position. The stream-ordered pool is what
        // makes that affordable -- plain cudaMalloc is tens of microseconds a
        // call, and a wide column needs tens of them. Widening later appends
        // to `bases` without disturbing anything already allocated.
        const std::size_t bytes = rows_ * sizeof(limb_t);
        c.bases.resize(c.nlimbs);
        for (std::size_t k = 0; k < c.nlimbs; ++k) {
            cuda(MallocAsync(&c.bases[k], bytes, 0));
            cuda(MemsetAsync(c.bases[k], 0, bytes, 0));
        }
        cuda(Malloc(&c.dev_bases, c.nlimbs * sizeof(limb_t *)));
        cuda(Memcpy(c.dev_bases, c.bases.data(), c.nlimbs * sizeof(limb_t *),
                    cudaMemcpyHostToDevice));
    }
    descriptors_stale_ = true;
}

// Shared tail of both reserve_for entry points: turn per-column exponent
// extents into reservations. Mirrors ColumnBlockMatrix::reserve_for, with one
// deliberate difference -- a column with nothing in it is still reserved,
// since the device cannot grow one later.
void
CudaColumnBlockMatrix::reserve_from_extents(const long long *low,
                                            const long long *high,
                                            const char *any, std::size_t count)
{
    const std::size_t headroom =
        ceil_log2(count < 1 ? 1 : count) + 1;  // +1 for the sign
    for (std::size_t j = 0; j < cols_; ++j) {
        if (!any[j]) {
            reserve_column(j, 0, limbs::kLimbBits);
            continue;
        }
        reserve_column(j, static_cast<int>(low[j]),
                       static_cast<std::size_t>(high[j] - low[j]) + headroom);
    }
}

void
CudaColumnBlockMatrix::reserve_for(const double *b, std::size_t count,
                                   std::size_t col_stride)
{
    if (0 == cols_) return;
    const std::size_t stride = col_stride ? col_stride : rows_;
    std::vector<long long> low(cols_, 0), high(cols_, 0);
    std::vector<char> any(cols_, 0);

    for (std::size_t j = 0; j < cols_; ++j) {
        // A column with no scale yet needs the exact true-ulp minimum, so the
        // floor is one nothing can clear.
        const kernels::Survey sv = kernels::survey()(
            b + j * stride, rows_, std::numeric_limits<long long>::max());
        if (sv.nonfinite) {
            throw std::domain_error(
                "cbfp: cannot reserve from a non-finite value");
        }
        if (!sv.any) continue;
        low[j] = sv.min_exponent;
        high[j] = sv.max_top;
        any[j] = 1;
    }
    reserve_from_extents(low.data(), high.data(), any.data(), count);
}

void
CudaColumnBlockMatrix::reserve_for_device(const double *b, std::size_t count,
                                          std::size_t col_stride)
{
    if (0 == cols_ || 0 == rows_) return;
    const std::size_t stride = col_stride ? col_stride : rows_;

    std::vector<DeviceSurvey> surveys(cols_);
    for (auto &s : surveys) s = DeviceSurvey{INT_MAX, INT_MIN, 0, 0};
    cuda(MemcpyAsync(survey_out_, surveys.data(), cols_ * sizeof(DeviceSurvey),
                     cudaMemcpyHostToDevice, st_compute_));
    survey_kernel<<<launch_grid(rows_, cols_), kBlock, 0, st_compute_>>>(
        b, rows_, stride, static_cast<DeviceSurvey *>(survey_out_));
    cuda(MemcpyAsync(surveys.data(), survey_out_, cols_ * sizeof(DeviceSurvey),
                     cudaMemcpyDeviceToHost, st_compute_));
    synchronize();

    std::vector<long long> low(cols_, 0), high(cols_, 0);
    std::vector<char> any(cols_, 0);
    for (std::size_t j = 0; j < cols_; ++j) {
        if (surveys[j].nonfinite) {
            throw std::domain_error(
                "cbfp: cannot reserve from a non-finite value");
        }
        if (!surveys[j].any) continue;
        low[j] = surveys[j].min_exponent;
        high[j] = surveys[j].max_top;
        any[j] = 1;
    }
    reserve_from_extents(low.data(), high.data(), any.data(), count);
}

void
CudaColumnBlockMatrix::reserve_like(const ColumnBlockMatrix &cpu)
{
    if (cpu.rows() != rows_ || cpu.cols() != cols_) {
        throw std::runtime_error(
            "cbfp: reserve_like requires matching dimensions");
    }
    for (std::size_t j = 0; j < cols_; ++j) {
        reserve_column(j, cpu.column_exponent(j), cpu.column_bit_width(j));
    }
}

int
CudaColumnBlockMatrix::column_occupancy(std::size_t j) const
{
    if (j >= cols_) throw std::out_of_range("cbfp: column index out of range");
    synchronize();
    const_cast<CudaColumnBlockMatrix *>(this)->harvest_occupancy();
    return cols_state_[j].max_limb_used;
}

// The device staging accumulates across launches, so the mapped array already
// holds the high-water mark; this only copies it where the column keeps it.
void
CudaColumnBlockMatrix::harvest_occupancy()
{
    if (nullptr == occupancy_host_) return;
    for (std::size_t j = 0; j < cols_; ++j) {
        cols_state_[j].max_limb_used = occupancy_host_[j];
    }
}

std::size_t
CudaColumnBlockMatrix::memory_bytes() const
{
    std::size_t total = 0;
    for (const auto &c : cols_state_) {
        total += c.nlimbs * rows_ * sizeof(limb_t);
    }
    return total;
}

void
CudaColumnBlockMatrix::sync_descriptors()
{
    if (!descriptors_stale_) return;
    std::vector<ColumnDesc> host(cols_);
    for (std::size_t j = 0; j < cols_; ++j) {
        host[j].bases = cols_state_[j].dev_bases;
        host[j].exponent = cols_state_[j].exponent;
        host[j].nlimbs = static_cast<unsigned>(cols_state_[j].nlimbs);
    }
    cuda(Memcpy(descriptors_, host.data(), cols_ * sizeof(ColumnDesc),
                cudaMemcpyHostToDevice));
    descriptors_stale_ = false;
}

void
CudaColumnBlockMatrix::ensure_slots(std::size_t words)
{
    if (slot_words_ >= words) return;
    synchronize();
    for (Slot &s : slots_) {
        if (s.pinned) cudaFreeHost(s.pinned);
        if (s.device) cudaFree(s.device);
        s.pinned = nullptr;
        s.device = nullptr;
        cuda(
            HostAlloc(&s.pinned, words * sizeof(double), cudaHostAllocDefault));
        cuda(Malloc(&s.device, words * sizeof(double)));
        if (!s.ev_done) {
            cuda(EventCreateWithFlags(&s.ev_done, cudaEventDisableTiming));
            cuda(EventCreateWithFlags(&s.ev_copied, cudaEventDisableTiming));
        }
        s.in_flight = false;
    }
    slot_words_ = words;
}

// Validates a packed column-major batch on the host, using the same survey the
// CPU accumulator uses -- which is the AVX-512 one where available, at ~0.12
// ns/elem. Cheap enough to hide entirely behind the transfer it runs against.
void
CudaColumnBlockMatrix::validate_host_survey(const double *packed)
{
    for (std::size_t j = 0; j < cols_; ++j) {
        const Column &c = cols_state_[j];
        // Passing the column's own exponent as the floor lets the survey stop
        // at a lower bound whenever no rescale could be due, and only pay for
        // the exact true-ulp minimum when it might be.
        const kernels::Survey sc =
            kernels::survey()(packed + j * rows_, rows_, c.exponent);
        if (sc.nonfinite) {
            throw std::domain_error(
                "cbfp: cannot accumulate a non-finite value");
        }
        if (!sc.any) continue;
        if (sc.min_exponent < c.exponent) {
            std::ostringstream os;
            os << "cbfp: column " << j << " was reserved at exponent "
               << c.exponent << " but these values need " << sc.min_exponent
               << "; pre-size it with reserve_column or reserve_like";
            throw std::runtime_error(os.str());
        }
        const long long width = sc.max_top - c.exponent;
        if (width > static_cast<long long>(c.nlimbs * limbs::kLimbBits)) {
            std::ostringstream os;
            os << "cbfp: column " << j << " was reserved at "
               << c.nlimbs * limbs::kLimbBits << " bits but these values reach "
               << width;
            throw std::runtime_error(os.str());
        }
    }
}

// The host path runs the transfer and the survey against each other: stage into
// pinned memory, start the copy, then survey that same buffer while the DMA is
// in flight. The survey is pure host work on host memory, so it costs nothing
// the transfer was not already going to spend.
//
// Validation therefore still happens *before* the accumulate is launched, so a
// bad batch is rejected without having touched the accumulator -- which a
// device-side survey cannot do without a round-trip that drains the pipeline.
void
CudaColumnBlockMatrix::add_matrix_col_major(const double *b,
                                            std::size_t col_stride)
{
    if (0 == rows_ || 0 == cols_) return;
    const std::size_t stride = col_stride ? col_stride : rows_;
    ensure_slots(rows_ * cols_);

    Slot &s = slots_[slot_];
    // A slot cannot be refilled until the kernel that last read it is done.
    if (s.in_flight) {
        cuda(EventSynchronize(s.ev_done));
        s.in_flight = false;
    }

    for (std::size_t j = 0; j < cols_; ++j) {
        std::memcpy(s.pinned + j * rows_, b + j * stride,
                    rows_ * sizeof(double));
    }
    cuda(MemcpyAsync(s.device, s.pinned, rows_ * cols_ * sizeof(double),
                     cudaMemcpyHostToDevice, st_copy_));
    cuda(EventRecord(s.ev_copied, st_copy_));

    // Runs against the copy above, not after it.
    validate_host_survey(s.pinned);

    // The accumulate waits on this slot's copy, but nothing else does, so the
    // next batch's transfer proceeds on the copy engine meanwhile.
    cuda(StreamWaitEvent(st_compute_, s.ev_copied, 0));
    launch_accumulate(s.device, rows_);
    cuda(EventRecord(s.ev_done, st_compute_));
    s.in_flight = true;
    slot_ ^= 1;
}

void
CudaColumnBlockMatrix::add_matrix_col_major_device(const double *b,
                                                   std::size_t col_stride)
{
    if (0 == rows_ || 0 == cols_) return;
    accumulate_device(b, col_stride ? col_stride : rows_);
}

void
CudaColumnBlockMatrix::synchronize() const
{
    cuda(StreamSynchronize(st_copy_));
    cuda(StreamSynchronize(st_compute_));
}

void
CudaColumnBlockMatrix::accumulate_device(const double *b,
                                         std::size_t col_stride)
{
    for (std::size_t j = 0; j < cols_; ++j) {
        if (!cols_state_[j].reserved) {
            throw std::runtime_error(
                "cbfp: every device column must be reserved before "
                "accumulation");
        }
    }
    sync_descriptors();

    const dim3 grid = launch_grid(rows_, cols_);

    // First pass: learn each column's exponent range, and reject anything the
    // reservation cannot hold. The decisions the CPU makes by rescaling and
    // widening are errors here, because neither is possible mid-launch.
    std::vector<DeviceSurvey> surveys(cols_);
    for (auto &s : surveys) s = DeviceSurvey{INT_MAX, INT_MIN, 0, 0};
    cuda(MemcpyAsync(survey_out_, surveys.data(), cols_ * sizeof(DeviceSurvey),
                     cudaMemcpyHostToDevice, st_compute_));

    survey_kernel<<<grid, kBlock, 0, st_compute_>>>(
        b, rows_, col_stride, static_cast<DeviceSurvey *>(survey_out_));
    // This is the drain the host path avoids: with the input already on the
    // device there is nothing to survey on the host, so the verdict has to come
    // back before the accumulate can be allowed to run.
    cuda(MemcpyAsync(surveys.data(), survey_out_, cols_ * sizeof(DeviceSurvey),
                     cudaMemcpyDeviceToHost, st_compute_));
    synchronize();

    for (std::size_t j = 0; j < cols_; ++j) {
        if (surveys[j].nonfinite) {
            throw std::domain_error(
                "cbfp: cannot accumulate a non-finite value");
        }
        if (!surveys[j].any) continue;
        const Column &c = cols_state_[j];
        if (surveys[j].min_exponent < c.exponent) {
            std::ostringstream os;
            os << "cbfp: column " << j << " was reserved at exponent "
               << c.exponent << " but these values need "
               << surveys[j].min_exponent
               << "; pre-size it with reserve_column or reserve_like";
            throw std::runtime_error(os.str());
        }
        const long long width = surveys[j].max_top - c.exponent;
        if (width > static_cast<long long>(c.nlimbs * limbs::kLimbBits)) {
            std::ostringstream os;
            os << "cbfp: column " << j << " was reserved at "
               << c.nlimbs * limbs::kLimbBits << " bits but these values reach "
               << width;
            throw std::runtime_error(os.str());
        }
    }

    launch_accumulate(b, col_stride);
}

// Queues the accumulate. Asynchronous: the caller is not blocked, so batches
// pipeline against each other without the caller doing anything.
void
CudaColumnBlockMatrix::launch_accumulate(const double *b,
                                         std::size_t col_stride)
{
    for (std::size_t j = 0; j < cols_; ++j) {
        if (!cols_state_[j].reserved) {
            throw std::runtime_error(
                "cbfp: every device column must be reserved before "
                "accumulation");
        }
    }
    sync_descriptors();

    const dim3 grid = launch_grid(rows_, cols_);

    accumulate_kernel<64><<<grid, kBlock, 0, st_compute_>>>(
        static_cast<const ColumnDesc *>(descriptors_), b, rows_, col_stride,
        occupancy_device_, occupancy_host_, ticket_,
        static_cast<unsigned>(cols_));
}

std::vector<limb_t>
CudaColumnBlockMatrix::entry_limbs(std::size_t i, std::size_t j) const
{
    check_index(i, j);
    synchronize();
    const Column &c = cols_state_[j];
    std::vector<limb_t> v(c.nlimbs);
    // One small copy per limb position. This is a readback path, not a hot
    // one; a bulk download would gather whole limb arrays instead.
    for (std::size_t k = 0; k < c.nlimbs; ++k) {
        cuda(Memcpy(&v[k], c.bases[k] + i, sizeof(limb_t),
                    cudaMemcpyDeviceToHost));
    }
    return v;
}

std::vector<limb_t>
CudaColumnBlockMatrix::download_column(std::size_t j) const
{
    if (j >= cols_) throw std::out_of_range("cbfp: column index out of range");
    synchronize();
    const Column &c = cols_state_[j];
    std::vector<limb_t> out(c.nlimbs * rows_);
    for (std::size_t k = 0; k < c.nlimbs; ++k) {
        cuda(Memcpy(out.data() + k * rows_, c.bases[k], rows_ * sizeof(limb_t),
                    cudaMemcpyDeviceToHost));
    }
    return out;
}

}  // namespace cbfp

#undef cuda
