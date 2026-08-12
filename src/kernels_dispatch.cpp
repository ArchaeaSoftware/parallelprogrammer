// Compiled at baseline. Selects a kernel variant once, from the running CPU.
#include <cstdlib>
#include <cstring>

#include "kernels.hpp"

namespace cbfp {
namespace kernels {
namespace {

AccumulateFn select_accumulate()
{
    // Escape hatch so both variants can be measured and cross-checked from a
    // single binary.
    if (const char* forced = std::getenv("CBFP_KERNEL")) {
        if (std::strcmp(forced, "scalar") == 0) return accumulate_scalar;
    }
#if defined(CBFP_HAVE_AVX512)
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx512f") &&
        __builtin_cpu_supports("avx512dq") &&
        __builtin_cpu_supports("avx512vl")) {
        return accumulate_avx512;
    }
#endif
    return accumulate_scalar;
}

}  // namespace

AccumulateFn accumulate()
{
    static const AccumulateFn fn = select_accumulate();
    return fn;
}

const char* accumulate_name()
{
    return accumulate() == accumulate_scalar ? "scalar" : "avx512";
}

}  // namespace kernels
}  // namespace cbfp
