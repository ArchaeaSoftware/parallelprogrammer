#include <cuda_runtime.h>

#include <climits>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "cbfp/column_accumulator.hpp"
#include "cbfp/cuda_accumulator.hpp"
#include "kernels.hpp"

namespace cbfp {
namespace {

using limbs::limb_t;

[[noreturn]] void cuda_fail(cudaError_t status, const char* call,
                            const char* file, int line)
{
    std::ostringstream os;
    os << "cbfp: " << call << " failed at " << file << ":" << line << ": "
       << cudaGetErrorName(status) << " -- " << cudaGetErrorString(status);
    throw std::runtime_error(os.str());
}

inline void cuda_check(cudaError_t status, const char* call, const char* file,
                       int line)
{
    if (status != cudaSuccess) cuda_fail(status, call, file, line);
}

// Every checked runtime call reads exactly like the call it makes, with one
// extra parenthesis: cuda(Malloc(&p, bytes)) is cudaMalloc(&p, bytes). The
// macro pastes the prefix back on for the call and stringizes it for the
// diagnostic, so a failure names the function rather than an opaque
// expression. A launch is checked with cuda(GetLastError()) and a blocking
// wait with cuda(DeviceSynchronize()), which paste the same way, so there is
// one mechanism rather than three.
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
constexpr unsigned kMaxGridX = 1024;

// What the first pass learns about a column, mirroring kernels::Scan. Reduced
// across the whole column with atomics, so every field is an atomic-friendly
// type rather than the host struct's bools.
struct DeviceScan {
    long long min_exponent;
    long long max_top;
    int any;
    int nonfinite;
};

struct ColumnDesc {
    // One base per limb position, so widening appends rather than moves.
    // A warp's threads all work different rows of the same limb column, so
    // this pointer load is a broadcast of one value, L1-resident for the life
    // of the kernel: the array is nlimbs * 8 bytes, 128 for a 16-limb column.
    limb_t* const* bases;
    long long exponent;
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
__device__ inline void split_device(double v, unsigned long long& mantissa,
                                    long long& exponent, long long& top,
                                    bool& negative, bool& nonfinite)
{
    const unsigned long long bits =
        static_cast<unsigned long long>(__double_as_longlong(v));
    const unsigned long long biased = (bits >> 52) & 0x7FFull;
    const unsigned long long frac = bits & ((1ull << 52) - 1);

    unsigned long long m = frac;
    if (biased != 0) m |= 1ull << 52;
    long long e = static_cast<long long>(biased < 1 ? 1 : biased) - 1075;

    negative = (bits >> 63) != 0;
    nonfinite = (biased == 0x7FFull);
    // Normalizing to odd raises the exponent by the trailing zero count and
    // lowers the significand's width by the same amount, so `top` is
    // unaffected and never needs that count.
    top = e + (m == 0 ? 0 : 64 - __clzll(static_cast<long long>(m)));
    if (m != 0) {
        const int tz = __ffsll(static_cast<long long>(m)) - 1;
        m >>= tz;
        e += tz;
    }
    mantissa = m;
    exponent = e;
}

// The CPU scan splits into two passes so the significand can be skipped once a
// column's scale is settled. That does not pay here: this pass is bound by
// reading the column, and the trailing-zero count is a few ALU ops on a value
// already in registers. One exact pass is both simpler and cheaper.
__global__ void scan_kernel(const double* __restrict__ values,
                            unsigned long long rows,
                            unsigned long long col_stride,
                            DeviceScan* __restrict__ out)
{
    __shared__ long long s_min[kBlock];
    __shared__ long long s_max[kBlock];
    __shared__ int s_any[kBlock];
    __shared__ int s_bad[kBlock];

    const unsigned j = blockIdx.y;
    const double* col = values + static_cast<std::size_t>(j) * col_stride;

    long long tmin = LLONG_MAX;
    long long tmax = LLONG_MIN;
    int tany = 0;
    int tbad = 0;

    const unsigned long long step =
        static_cast<unsigned long long>(gridDim.x) * blockDim.x;
    for (unsigned long long i =
             static_cast<unsigned long long>(blockIdx.x) * blockDim.x +
             threadIdx.x;
         i < rows; i += step) {
        unsigned long long m;
        long long e, top;
        bool neg, bad;
        split_device(col[i], m, e, top, neg, bad);
        if (bad) tbad = 1;
        if (m != 0) {
            tany = 1;
            if (e < tmin) tmin = e;
            if (top > tmax) tmax = top;
        }
    }

    const unsigned t = threadIdx.x;
    s_min[t] = tmin;
    s_max[t] = tmax;
    s_any[t] = tany;
    s_bad[t] = tbad;
    __syncthreads();
    for (unsigned s = blockDim.x / 2; s > 0; s >>= 1) {
        if (t < s) {
            if (s_min[t + s] < s_min[t]) s_min[t] = s_min[t + s];
            if (s_max[t + s] > s_max[t]) s_max[t] = s_max[t + s];
            s_any[t] |= s_any[t + s];
            s_bad[t] |= s_bad[t + s];
        }
        __syncthreads();
    }

    // One atomic per block rather than one per thread.
    if (t == 0) {
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
__global__ void accumulate_kernel(const ColumnDesc* __restrict__ cols,
                                  const double* __restrict__ values,
                                  unsigned long long rows,
                                  unsigned long long col_stride)
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
    const double* col = values + static_cast<std::size_t>(j) * col_stride;

    const unsigned long long step =
        static_cast<unsigned long long>(gridDim.x) * blockDim.x;
    for (unsigned long long i =
             static_cast<unsigned long long>(blockIdx.x) * blockDim.x +
             threadIdx.x;
         i < rows; i += step) {
        unsigned long long m;
        long long e, top;
        bool neg, bad;
        split_device(col[i], m, e, top, neg, bad);
        if (m == 0) continue;

        // The host has already checked that every value fits the column it was
        // reserved for, so the shift cannot be negative.
        const unsigned long long shift =
            static_cast<unsigned long long>(e - c.exponent);
        const unsigned off = static_cast<unsigned>(shift / kRadix);
        const unsigned bit = static_cast<unsigned>(shift % kRadix);

        // A 53-bit significand at intra-limb offset `bit` spans two limbs at
        // radix 64 and at radix 52 alike.
        unsigned long long lo, hi;
        if (bit == 0) {
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
            limb_t* dst = c.bases[p] + i;
            const unsigned long long x = *dst;
            if (neg) {
                const unsigned long long d = x - a;
                const unsigned long long b1 = (x < a) ? 1ull : 0ull;
                const unsigned long long d2 = d - carry;
                const unsigned long long b2 = (d < carry) ? 1ull : 0ull;
                *dst = d2;
                carry = b1 | b2;
            } else {
                const unsigned long long s = x + a;
                const unsigned long long c1 = (s < x) ? 1ull : 0ull;
                const unsigned long long s2 = s + carry;
                const unsigned long long c2 = (s2 < s) ? 1ull : 0ull;
                *dst = s2;
                carry = c1 | c2;
            }
            if (carry == 0 && p >= off + 1) break;
        }
    }
}

std::size_t limbs_for_bits(std::size_t bits)
{
    return (bits + limbs::kLimbBits - 1) / limbs::kLimbBits;
}

}  // namespace

bool cuda_available()
{
    int n = 0;
    return cudaGetDeviceCount(&n) == cudaSuccess && n > 0;
}

CudaColumnBlockMatrix::CudaColumnBlockMatrix(std::size_t rows, std::size_t cols)
    : rows_(rows), cols_(cols), cols_state_(cols)
{
    if (!cuda_available()) {
        throw std::runtime_error("cbfp: no usable CUDA device");
    }
    cudaStream_t st;
    cuda(StreamCreate(&st));
    stream_ = st;

    if (cols_ != 0) {
        cuda(Malloc(&descriptors_, cols_ * sizeof(ColumnDesc)));
        cuda(Malloc(&scan_out_, cols_ * sizeof(DeviceScan)));
    }
}

CudaColumnBlockMatrix::~CudaColumnBlockMatrix()
{
    // Deliberately unchecked: a destructor must not throw, and there is
    // nothing useful to do about a failed free during teardown.
    for (auto& c : cols_state_) {
        for (limb_t* p : c.bases) cudaFreeAsync(p, 0);
        cudaFree(c.dev_bases);
    }
    cudaFree(descriptors_);
    cudaFree(scan_out_);
    for (Slot& s : slots_) {
        if (s.pinned) cudaFreeHost(s.pinned);
        if (s.device) cudaFree(s.device);
        if (s.done) cudaEventDestroy(static_cast<cudaEvent_t>(s.done));
    }
    if (stream_) cudaStreamDestroy(static_cast<cudaStream_t>(stream_));
}

void CudaColumnBlockMatrix::check_index(std::size_t i, std::size_t j) const
{
    if (i >= rows_ || j >= cols_) {
        throw std::out_of_range("cbfp: matrix index out of range");
    }
}

int CudaColumnBlockMatrix::column_exponent(std::size_t j) const
{
    if (j >= cols_) throw std::out_of_range("cbfp: column index out of range");
    return cols_state_[j].exponent;
}

std::size_t CudaColumnBlockMatrix::column_limbs(std::size_t j) const
{
    if (j >= cols_) throw std::out_of_range("cbfp: column index out of range");
    return cols_state_[j].nlimbs;
}

void CudaColumnBlockMatrix::reserve_column(std::size_t j, int exponent,
                                           std::size_t bits)
{
    if (j >= cols_) throw std::out_of_range("cbfp: column index out of range");
    Column& c = cols_state_[j];
    if (c.reserved) {
        throw std::runtime_error(
            "cbfp: a device column may only be reserved once");
    }

    c.exponent = exponent;
    c.nlimbs = limbs_for_bits(bits) < 1 ? 1 : limbs_for_bits(bits);
    c.reserved = true;

    if (rows_ != 0) {
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
        cuda(Malloc(&c.dev_bases, c.nlimbs * sizeof(limb_t*)));
        cuda(Memcpy(c.dev_bases, c.bases.data(), c.nlimbs * sizeof(limb_t*),
                    cudaMemcpyHostToDevice));
    }
    descriptors_stale_ = true;
}

void CudaColumnBlockMatrix::reserve_like(const ColumnBlockMatrix& cpu)
{
    if (cpu.rows() != rows_ || cpu.cols() != cols_) {
        throw std::runtime_error(
            "cbfp: reserve_like requires matching dimensions");
    }
    for (std::size_t j = 0; j < cols_; ++j) {
        reserve_column(j, cpu.column_exponent(j), cpu.column_bit_width(j));
    }
}

std::size_t CudaColumnBlockMatrix::memory_bytes() const
{
    std::size_t total = 0;
    for (const auto& c : cols_state_) {
        total += c.nlimbs * rows_ * sizeof(limb_t);
    }
    return total;
}

void CudaColumnBlockMatrix::sync_descriptors()
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

void CudaColumnBlockMatrix::ensure_slots(std::size_t words)
{
    if (slot_words_ >= words) return;
    synchronize();
    for (Slot& s : slots_) {
        if (s.pinned) cudaFreeHost(s.pinned);
        if (s.device) cudaFree(s.device);
        s.pinned = nullptr;
        s.device = nullptr;
        cuda(
            HostAlloc(&s.pinned, words * sizeof(double), cudaHostAllocDefault));
        cuda(Malloc(&s.device, words * sizeof(double)));
        if (!s.done) {
            cudaEvent_t e;
            cuda(EventCreateWithFlags(&e, cudaEventDisableTiming));
            s.done = e;
        }
        s.in_flight = false;
    }
    slot_words_ = words;
}

// Validates a packed column-major batch on the host, using the same scan the
// CPU accumulator uses -- which is the AVX-512 one where available, at ~0.12
// ns/elem. Cheap enough to hide entirely behind the transfer it runs against.
void CudaColumnBlockMatrix::validate_host_scan(const double* packed)
{
    for (std::size_t j = 0; j < cols_; ++j) {
        const Column& c = cols_state_[j];
        // Passing the column's own exponent as the floor lets the scan stop at
        // a lower bound whenever no rescale could be due, and only pay for the
        // exact true-ulp minimum when it might be.
        const kernels::Scan sc =
            kernels::scan()(packed + j * rows_, rows_, c.exponent);
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

// The host path runs the transfer and the scan against each other: stage into
// pinned memory, start the copy, then scan that same buffer while the DMA is
// in flight. The scan is pure host work on host memory, so it costs nothing
// the transfer was not already going to spend.
//
// Validation therefore still happens *before* the accumulate is launched, so a
// bad batch is rejected without having touched the accumulator -- which a
// device-side scan cannot do without a round-trip that drains the pipeline.
void CudaColumnBlockMatrix::add_matrix_col_major(const double* b,
                                                 std::size_t col_stride)
{
    if (rows_ == 0 || cols_ == 0) return;
    const std::size_t stride = col_stride ? col_stride : rows_;
    ensure_slots(rows_ * cols_);

    Slot& s = slots_[slot_];
    // A slot cannot be refilled until the kernel that last read it is done.
    if (s.in_flight) {
        cuda(EventSynchronize(static_cast<cudaEvent_t>(s.done)));
        s.in_flight = false;
    }

    for (std::size_t j = 0; j < cols_; ++j) {
        std::memcpy(s.pinned + j * rows_, b + j * stride,
                    rows_ * sizeof(double));
    }
    cuda(MemcpyAsync(s.device, s.pinned, rows_ * cols_ * sizeof(double),
                     cudaMemcpyHostToDevice,
                     static_cast<cudaStream_t>(stream_)));

    // Runs against the copy above, not after it.
    validate_host_scan(s.pinned);

    launch_accumulate(s.device, rows_);
    cuda(EventRecord(static_cast<cudaEvent_t>(s.done),
                     static_cast<cudaStream_t>(stream_)));
    s.in_flight = true;
    slot_ ^= 1;
}

void CudaColumnBlockMatrix::add_matrix_col_major_device(const double* b,
                                                        std::size_t col_stride)
{
    if (rows_ == 0 || cols_ == 0) return;
    accumulate_device(b, col_stride ? col_stride : rows_);
}

void CudaColumnBlockMatrix::synchronize() const
{
    cuda(StreamSynchronize(static_cast<cudaStream_t>(stream_)));
}

void CudaColumnBlockMatrix::accumulate_device(const double* b,
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

    const unsigned gx =
        static_cast<unsigned>((rows_ + kBlock - 1) / kBlock > kMaxGridX
                                  ? kMaxGridX
                                  : (rows_ + kBlock - 1) / kBlock);
    const dim3 grid(gx == 0 ? 1 : gx, static_cast<unsigned>(cols_));

    // First pass: learn each column's exponent range, and reject anything the
    // reservation cannot hold. The decisions the CPU makes by rescaling and
    // widening are errors here, because neither is possible mid-launch.
    std::vector<DeviceScan> scans(cols_);
    for (auto& s : scans) s = DeviceScan{LLONG_MAX, LLONG_MIN, 0, 0};
    cuda(MemcpyAsync(scan_out_, scans.data(), cols_ * sizeof(DeviceScan),
                     cudaMemcpyHostToDevice,
                     static_cast<cudaStream_t>(stream_)));

    scan_kernel<<<grid, kBlock, 0, static_cast<cudaStream_t>(stream_)>>>(
        b, rows_, col_stride, static_cast<DeviceScan*>(scan_out_));
    cuda(GetLastError());
    // This is the drain the host path avoids: with the input already on the
    // device there is nothing to scan on the host, so the verdict has to come
    // back before the accumulate can be allowed to run.
    cuda(MemcpyAsync(scans.data(), scan_out_, cols_ * sizeof(DeviceScan),
                     cudaMemcpyDeviceToHost,
                     static_cast<cudaStream_t>(stream_)));
    synchronize();

    for (std::size_t j = 0; j < cols_; ++j) {
        if (scans[j].nonfinite) {
            throw std::domain_error(
                "cbfp: cannot accumulate a non-finite value");
        }
        if (!scans[j].any) continue;
        const Column& c = cols_state_[j];
        if (scans[j].min_exponent < c.exponent) {
            std::ostringstream os;
            os << "cbfp: column " << j << " was reserved at exponent "
               << c.exponent << " but these values need "
               << scans[j].min_exponent
               << "; pre-size it with reserve_column or reserve_like";
            throw std::runtime_error(os.str());
        }
        const long long width = scans[j].max_top - c.exponent;
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
void CudaColumnBlockMatrix::launch_accumulate(const double* b,
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

    const unsigned gx =
        static_cast<unsigned>((rows_ + kBlock - 1) / kBlock > kMaxGridX
                                  ? kMaxGridX
                                  : (rows_ + kBlock - 1) / kBlock);
    const dim3 grid(gx == 0 ? 1 : gx, static_cast<unsigned>(cols_));

    accumulate_kernel<64>
        <<<grid, kBlock, 0, static_cast<cudaStream_t>(stream_)>>>(
            static_cast<const ColumnDesc*>(descriptors_), b, rows_, col_stride);
    cuda(GetLastError());
}

std::vector<limb_t> CudaColumnBlockMatrix::entry_limbs(std::size_t i,
                                                       std::size_t j) const
{
    check_index(i, j);
    synchronize();
    const Column& c = cols_state_[j];
    std::vector<limb_t> v(c.nlimbs);
    // One small copy per limb position. This is a readback path, not a hot
    // one; a bulk download would gather whole limb arrays instead.
    for (std::size_t k = 0; k < c.nlimbs; ++k) {
        cuda(Memcpy(&v[k], c.bases[k] + i, sizeof(limb_t),
                    cudaMemcpyDeviceToHost));
    }
    return v;
}

std::vector<limb_t> CudaColumnBlockMatrix::download_column(std::size_t j) const
{
    if (j >= cols_) throw std::out_of_range("cbfp: column index out of range");
    synchronize();
    const Column& c = cols_state_[j];
    std::vector<limb_t> out(c.nlimbs * rows_);
    for (std::size_t k = 0; k < c.nlimbs; ++k) {
        cuda(Memcpy(out.data() + k * rows_, c.bases[k], rows_ * sizeof(limb_t),
                    cudaMemcpyDeviceToHost));
    }
    return out;
}

}  // namespace cbfp

#undef cuda
