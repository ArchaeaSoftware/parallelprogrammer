// Elementwise AVX-512 operators for the fusion experiment.
//
// Every operator holds its constants as __m512 members initialised from
// runtime values, so the compiler cannot constant-fold the arithmetic but
// can still hoist the broadcasts out of the loop. Each op is a distinct
// type, which is what makes the variadic composition inlinable.
#pragma once
#include <immintrin.h>
#include <cstddef>

// max(v, 0) -- cheapest possible link: one constant, one instruction.
struct Relu {
    __m512 z;
    explicit Relu(float zero) : z(_mm512_set1_ps(zero)) {}
    __m512 operator()(__m512 v) const { return _mm512_max_ps(v, z); }
    static const char* name() { return "relu"; }
};

// v * k -- baseline multiply: one constant, one instruction.
struct Scale {
    __m512 k;
    explicit Scale(float k_) : k(_mm512_set1_ps(k_)) {}
    __m512 operator()(__m512 v) const { return _mm512_mul_ps(v, k); }
    static const char* name() { return "scale"; }
};

// v * a + b -- same instruction count as Scale, twice the register cost.
struct FmaBias {
    __m512 a, b;
    FmaBias(float a_, float b_) : a(_mm512_set1_ps(a_)), b(_mm512_set1_ps(b_)) {}
    __m512 operator()(__m512 v) const { return _mm512_fmadd_ps(v, a, b); }
    static const char* name() { return "fma_bias"; }
};

// min(max(v, lo), hi) -- two constants, two instructions.
struct Clamp {
    __m512 lo, hi;
    Clamp(float lo_, float hi_) : lo(_mm512_set1_ps(lo_)), hi(_mm512_set1_ps(hi_)) {}
    __m512 operator()(__m512 v) const {
        return _mm512_min_ps(_mm512_max_ps(v, lo), hi);
    }
    static const char* name() { return "clamp"; }
};

// Horner-form polynomial of degree D: D+1 constants, D dependent FMAs.
// The sweep knob -- raising D lengthens the serial chain and consumes one
// more coefficient register per degree.
template <int D>
struct Poly {
    __m512 c[D + 1];
    explicit Poly(const float* coeffs) {
        for (int i = 0; i <= D; ++i) c[i] = _mm512_set1_ps(coeffs[i]);
    }
    __m512 operator()(__m512 v) const {
        __m512 r = c[D];
        for (int i = D - 1; i >= 0; --i) r = _mm512_fmadd_ps(r, v, c[i]);
        return r;
    }
    static const char* name() { return "poly"; }
};

// 16-entry table lookup via vpermps -- contends for the shuffle port rather
// than the FMA ports. Shape of the NF4 dequantisation codebook lookup.
struct PermuteLut {
    __m512 table;
    __m512i mask;
    explicit PermuteLut(const float* tbl)
        : table(_mm512_loadu_ps(tbl)), mask(_mm512_set1_epi32(15)) {}
    __m512 operator()(__m512 v) const {
        __m512i idx = _mm512_and_si512(_mm512_cvttps_epi32(v), mask);
        return _mm512_permutexvar_ps(idx, table);
    }
    static const char* name() { return "permute_lut"; }
};

// rsqrt with one Newton-Raphson refinement -- long dependent chain,
// few constants. Requires positive input.
struct RsqrtNR {
    __m512 half, three;
    RsqrtNR(float h, float t) : half(_mm512_set1_ps(h)), three(_mm512_set1_ps(t)) {}
    __m512 operator()(__m512 v) const {
        __m512 x  = _mm512_rsqrt14_ps(v);
        __m512 vx = _mm512_mul_ps(v, x);
        __m512 t  = _mm512_fnmadd_ps(vx, x, three);   // 3 - v*x*x
        return _mm512_mul_ps(_mm512_mul_ps(x, half), t);
    }
    static const char* name() { return "rsqrt_nr"; }
};
