// Compiled at baseline. Selects a kernel variant once, from the running CPU.
#include <cstdlib>
#include <cstring>

#include "kernels.hpp"

namespace cbfp {
namespace kernels {
namespace {

bool use_avx512()
{
    // Escape hatch so both variants can be measured and cross-checked from a
    // single binary.
    if (const char* forced = std::getenv("CBFP_KERNEL")) {
        if (std::strcmp(forced, "scalar") == 0) return false;
    }
#if defined(CBFP_HAVE_AVX512)
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

AccumulateFn accumulate()
{
#if defined(CBFP_HAVE_AVX512)
    static const AccumulateFn fn =
        use_avx512() ? accumulate_avx512 : accumulate_scalar;
#else
    static const AccumulateFn fn = accumulate_scalar;
#endif
    return fn;
}

ScanFn scan()
{
#if defined(CBFP_HAVE_AVX512)
    static const ScanFn fn =
        use_avx512() ? scan_column_avx512 : scan_column_scalar;
#else
    static const ScanFn fn = scan_column_scalar;
#endif
    return fn;
}

const char* accumulate_name()
{
    return accumulate() == accumulate_scalar ? "scalar" : "avx512";
}

}  // namespace kernels
}  // namespace cbfp
