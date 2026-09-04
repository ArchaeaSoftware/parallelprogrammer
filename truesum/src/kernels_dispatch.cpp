// Compiled at baseline. Selects a kernel variant once, from the running CPU.
#include <cstdlib>
#include <cstring>

#include "kernels.hpp"

namespace truesum {
namespace kernels {
namespace {

bool
use_avx512()
{
    // Escape hatch so both variants can be measured and cross-checked from a
    // single binary.
    if (const char *forced = std::getenv("TRUESUM_KERNEL")) {
        if (0 == std::strcmp(forced, "scalar")) return false;
    }
#if defined(TRUESUM_HAVE_AVX512)
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx512f") &&
           __builtin_cpu_supports("avx512dq") &&
           __builtin_cpu_supports("avx512vl") &&
           __builtin_cpu_supports("avx512cd") &&
           // _mm512_popcnt_epi64, used for the vector trailing-zero count.
           __builtin_cpu_supports("avx512vpopcntdq");
#else
    return false;
#endif
}

}  // namespace

AccumulateFn
accumulate()
{
#if defined(TRUESUM_HAVE_AVX512)
    static const AccumulateFn fn =
        use_avx512() ? accumulate_avx512 : accumulate_scalar;
#else
    static const AccumulateFn fn = accumulate_scalar;
#endif
    return fn;
}

// The vector fold is what makes folding pay at all. Scalar, the fold moved
// less memory but gave up eight lanes to do it, and lost everywhere -- 0.19
// against 0.57 Gelem/s single-threaded, a flat 3x deficit at every shape. The
// traffic argument was sound; it just could not cover the cost of leaving
// AVX-512 behind.
AccumulateFoldFn
accumulate_fold()
{
#if defined(TRUESUM_HAVE_AVX512)
    static const AccumulateFoldFn fn =
        use_avx512() ? accumulate_fold_avx512 : accumulate_fold_scalar;
#else
    static const AccumulateFoldFn fn = accumulate_fold_scalar;
#endif
    return fn;
}

SurveyFn
survey()
{
#if defined(TRUESUM_HAVE_AVX512)
    static const SurveyFn fn =
        use_avx512() ? survey_column_avx512 : survey_column_scalar;
#else
    static const SurveyFn fn = survey_column_scalar;
#endif
    return fn;
}

const char *
accumulate_name()
{
    return accumulate() == accumulate_scalar ? "scalar" : "avx512";
}

}  // namespace kernels
}  // namespace truesum
