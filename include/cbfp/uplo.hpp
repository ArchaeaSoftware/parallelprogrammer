// Which triangle of a symmetric matrix is the stored one, following LAPACK's
// uplo.
//
// Lower keeps entries with i >= j, so column j holds rows j..n-1; Upper keeps
// i <= j, so column j holds rows 0..j. Either way a column's stored rows stay
// contiguous, which is what lets the vector and coalescing properties of the
// limb-major layout survive being made triangular.
//
// Its own header so that both accumulators can name it without the CUDA one
// having to include the CPU class it otherwise only forward-declares.
#pragma once

namespace cbfp {

enum class Uplo { Lower, Upper };

}  // namespace cbfp
