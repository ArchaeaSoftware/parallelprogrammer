// Storage for one limb position across every row of a matrix column.
//
// Entry `i`'s `k`-th limb lives at `limbs[k][i]`, so a LimbColumn is a slice at
// fixed significance spanning all rows. Note this is a limb *of* a column, not
// a column of the matrix.
#pragma once

#include <cstddef>
#include <cstdint>

namespace cbfp {

class LimbColumn {
public:
    // Rows are padded to this multiple so a full 512-bit vector access at the
    // final row block stays in bounds without a masked tail.
    static constexpr std::size_t kRowBlock = 8;

    // `index` is the limb position. It only skews the base address: limb
    // arrays are otherwise identical.
    LimbColumn(std::size_t rows, std::size_t index);
    ~LimbColumn();

    LimbColumn(LimbColumn &&other) noexcept;
    LimbColumn &operator=(LimbColumn &&other) noexcept;
    LimbColumn(const LimbColumn &) = delete;
    LimbColumn &operator=(const LimbColumn &) = delete;

    std::uint64_t *data() noexcept { return data_; }
    const std::uint64_t *data() const noexcept { return data_; }

private:
    void *alloc_ = nullptr;
    std::uint64_t *data_ = nullptr;
};

// Rows rounded up to a whole vector block.
std::size_t
padded_rows(std::size_t rows);

}  // namespace cbfp
