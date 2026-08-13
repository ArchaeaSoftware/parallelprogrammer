#include "cbfp/limb_column.hpp"

#include <cstdlib>
#include <cstring>
#include <new>
#include <utility>

namespace cbfp {
namespace {

// Successive limb arrays are offset by a multiple of the cache line so they do
// not land in the same cache set. Without this, allocations of a round size
// come back mutually 4KB-congruent and throughput becomes a lottery drawn by
// the allocator -- measured at 2.8x between outcomes. See docs/simd-design.md.
constexpr std::size_t kSkewCycle = 8;
constexpr std::size_t kSkewStep = 64;
constexpr std::size_t kAlign = 64;

}  // namespace

std::size_t padded_rows(std::size_t rows)
{
    return (rows + LimbColumn::kRowBlock - 1) & ~(LimbColumn::kRowBlock - 1);
}

LimbColumn::LimbColumn(std::size_t rows, std::size_t index)
{
    const std::size_t skew = (index % kSkewCycle) * kSkewStep;
    const std::size_t payload = padded_rows(rows) * sizeof(std::uint64_t);
    // aligned_alloc requires a size that is a multiple of the alignment.
    const std::size_t bytes =
        (payload + kSkewCycle * kSkewStep + kAlign - 1) & ~(kAlign - 1);

    alloc_ = std::aligned_alloc(kAlign, bytes);
    if (nullptr == alloc_) throw std::bad_alloc();
    std::memset(alloc_, 0, bytes);
    data_ = reinterpret_cast<std::uint64_t*>(static_cast<char*>(alloc_) + skew);
}

LimbColumn::~LimbColumn()
{
    std::free(alloc_);
}

LimbColumn::LimbColumn(LimbColumn&& other) noexcept
    : alloc_(other.alloc_), data_(other.data_)
{
    other.alloc_ = nullptr;
    other.data_ = nullptr;
}

LimbColumn& LimbColumn::operator=(LimbColumn&& other) noexcept
{
    if (this != &other) {
        std::free(alloc_);
        alloc_ = other.alloc_;
        data_ = other.data_;
        other.alloc_ = nullptr;
        other.data_ = nullptr;
    }
    return *this;
}

}  // namespace cbfp
