// Ways to apply a chain of elementwise operators to a buffer.
//
//   separate_passes      one full pass over memory per operator
//   fused_vecmajor       one pass; each vector runs its whole chain in turn
//   fused_vecmajor_seq   as above, unroll carried by an index_sequence pack
//   fused_opmajor_byval  one pass; each operator applied to all U accumulators
//   fused_opmajor        as above, accumulators captured by reference
//
// The pack carries the operators, which are heterogeneous. The unroll factor
// U is homogeneous repetition, carried either by a loop with a compile-time
// bound or by an index_sequence pack -- which turns out to matter at -O2,
// where GCC does not unroll and only the pack guarantees constant indices.
//
// Every path ends in a masked tail, so n need not be a multiple of 16.
#pragma once
#include <immintrin.h>
#include <cstddef>
#include <utility>
#include <array>

// Value loaded into inactive lanes of a masked tail load. Anything numerically
// dull will do; zero is a poor choice because it feeds infinities to rsqrt.
inline __m512 safe_fill() { return _mm512_set1_ps(1.0f); }

inline __mmask16 tail_mask(std::size_t rem) {
    return static_cast<__mmask16>((1u << rem) - 1u);
}

// Processes the final rem < 16 elements, if any.
template <typename... Ops>
inline void masked_tail(const float* in, float* out, std::size_t i,
                        std::size_t n, const Ops&... ops) {
    if (i >= n) return;
    const __mmask16 k = tail_mask(n - i);
    __m512 v = _mm512_mask_loadu_ps(safe_fill(), k, in + i);
    ((v = ops(v)), ...);
    _mm512_mask_storeu_ps(out + i, k, v);
}

// ---------------------------------------------------------------- separate

template <typename Op>
inline void one_pass(const float* in, float* out, std::size_t n, Op op) {
    std::size_t i = 0;
    for (; i + 16 <= n; i += 16)
        _mm512_storeu_ps(out + i, op(_mm512_loadu_ps(in + i)));
    masked_tail(in, out, i, n, op);
}

// Applies each operator as its own full pass, ping-ponging between buffers.
template <typename... Ops>
void separate_passes(const float* in, float* out, float* tmp, std::size_t n,
                     Ops... ops) {
    const float* src = in;
    float* buf[2] = { tmp, out };
    int b = 0;
    auto step = [&](auto op) {
        one_pass(src, buf[b], n, op);
        src = buf[b];
        b ^= 1;
    };
    (step(ops), ...);
    if (src != out) {
        std::size_t i = 0;
        for (; i + 16 <= n; i += 16)
            _mm512_storeu_ps(out + i, _mm512_loadu_ps(src + i));
        if (i < n) {
            const __mmask16 k = tail_mask(n - i);
            _mm512_mask_storeu_ps(out + i, k, _mm512_maskz_loadu_ps(k, src + i));
        }
    }
}

// ------------------------------------------------------------- vector-major

template <std::size_t U, typename... Ops>
void fused_vecmajor(const float* in, float* out, std::size_t n, Ops... ops) {
    constexpr std::size_t W = 16, BLK = W * U;
    std::size_t i = 0;
    for (; i + BLK <= n; i += BLK) {
        __m512 v[U];
        for (std::size_t u = 0; u < U; ++u) v[u] = _mm512_loadu_ps(in + i + W * u);
        for (std::size_t u = 0; u < U; ++u) ((v[u] = ops(v[u])), ...);
        for (std::size_t u = 0; u < U; ++u) _mm512_storeu_ps(out + i + W * u, v[u]);
    }
    for (; i + W <= n; i += W) {
        __m512 v = _mm512_loadu_ps(in + i);
        ((v = ops(v)), ...);
        _mm512_storeu_ps(out + i, v);
    }
    masked_tail(in, out, i, n, ops...);
}

// ------------------------- vector-major, unroll expanded by a parameter pack

template <typename... Ops>
inline __m512 chain(__m512 v, const Ops&... ops) {
    ((v = ops(v)), ...);
    return v;
}

template <std::size_t... I, typename... Ops>
inline void block_seq(const float* in, float* out, std::size_t i,
                      std::index_sequence<I...>, const Ops&... ops) {
    __m512 v[sizeof...(I)];
    ((v[I] = _mm512_loadu_ps(in + i + 16 * I)), ...);
    ((v[I] = chain(v[I], ops...)), ...);
    ((_mm512_storeu_ps(out + i + 16 * I, v[I])), ...);
}

template <std::size_t U, typename... Ops>
void fused_vecmajor_seq(const float* in, float* out, std::size_t n, Ops... ops) {
    constexpr std::size_t W = 16, BLK = W * U;
    std::size_t i = 0;
    for (; i + BLK <= n; i += BLK)
        block_seq(in, out, i, std::make_index_sequence<U>{}, ops...);
    for (; i + W <= n; i += W) {
        __m512 v = _mm512_loadu_ps(in + i);
        ((v = ops(v)), ...);
        _mm512_storeu_ps(out + i, v);
    }
    masked_tail(in, out, i, n, ops...);
}

// ------------------------------------------- op-major, accumulators by value

template <typename Op, std::size_t... I>
inline std::array<__m512, sizeof...(I)>
apply_all_byval(std::array<__m512, sizeof...(I)> v, const Op& op,
                std::index_sequence<I...>) {
    return { op(v[I])... };
}

template <std::size_t U, typename... Ops>
void fused_opmajor_byval(const float* in, float* out, std::size_t n, Ops... ops) {
    constexpr std::size_t W = 16, BLK = W * U;
    constexpr auto seq = std::make_index_sequence<U>{};
    std::size_t i = 0;
    for (; i + BLK <= n; i += BLK) {
        std::array<__m512, U> v;
        for (std::size_t u = 0; u < U; ++u) v[u] = _mm512_loadu_ps(in + i + W * u);
        ((v = apply_all_byval(v, ops, seq)), ...);
        for (std::size_t u = 0; u < U; ++u) _mm512_storeu_ps(out + i + W * u, v[u]);
    }
    for (; i + W <= n; i += W) {
        __m512 vv = _mm512_loadu_ps(in + i);
        ((vv = ops(vv)), ...);
        _mm512_storeu_ps(out + i, vv);
    }
    masked_tail(in, out, i, n, ops...);
}

// ----------------------------------------------------------------- op-major

template <std::size_t U, typename... Ops>
void fused_opmajor(const float* in, float* out, std::size_t n, Ops... ops) {
    constexpr std::size_t W = 16, BLK = W * U;
    std::size_t i = 0;
    for (; i + BLK <= n; i += BLK) {
        __m512 v[U];
        for (std::size_t u = 0; u < U; ++u) v[u] = _mm512_loadu_ps(in + i + W * u);
        auto apply_all = [&](auto op) {
            for (std::size_t u = 0; u < U; ++u) v[u] = op(v[u]);
        };
        (apply_all(ops), ...);
        for (std::size_t u = 0; u < U; ++u) _mm512_storeu_ps(out + i + W * u, v[u]);
    }
    for (; i + W <= n; i += W) {
        __m512 v = _mm512_loadu_ps(in + i);
        ((v = ops(v)), ...);
        _mm512_storeu_ps(out + i, v);
    }
    masked_tail(in, out, i, n, ops...);
}
