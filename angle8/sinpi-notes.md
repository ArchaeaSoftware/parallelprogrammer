# angle8: sin without the table, and with a table in registers

Notes from the `angle8-sinpi` exploration (commit 5e65ee4, one commit on top
of main, September 2026). The branch was dropped as not worth carrying; this
records what was tried and what came of it so it can be rebuilt if it ever is.
The commit stays reachable by hash until the reflog expires.

## The setup

`angle8_sin` only ever sees 65 distinct arguments. After the quarter-wave fold
the angle is k * pi/128 with k in 0..64; bit 6 of the raw byte mirrors the
index and bit 7 flips the sign:

```c
unsigned i = u & 63u;
unsigned k = (u & 64u) ? 64u - i : i;   /* 0..64 */
float    s = f(k);
return (u & 128u) ? 0.0f - s : s;       /* 0.0f - s keeps sin(pi) at +0 */
```

Every formulation shared that fold, so all of them are exactly odd, cos is
exactly even and `sin(a + pi/2)` is bit-identical to `cos(a)`. They differed
only in how they got from k to a float.

## Three scalar formulations

| | how | constants | worst error, plain | worst error, fma |
|---|---|---|---|---|
| table | 65-entry table, correctly rounded (current main) | 65 floats | 0.5 ulp | 0.5 ulp |
| poly | degree-9 odd polynomial in x = k/128 | 5 floats | 1.59 ulp | 1.29 ulp |
| split | k = 8 hi + lo; sin(hi) cos(lo) + cos(hi) sin(lo) | 25 floats | 1.53 ulp | 0.97 ulp |

Scalar timings, Ryzen 7 7700X, gcc 9, `-O2`, dependency-chained loop as in
`bench.c`: table 0.55 ns, split 0.75 ns, poly 1.06 ns, libm `sinf` 2.6 ns.

### Polynomial

sinpi(x) = x * q(x^2), q of degree 4 in x^2, fitted to the 65 reachable points
(discrete minimax in relative error, then a coordinate-descent search over
nearby float coefficients scoring plain and fused evaluation together, with
sinpi(1/2) == 1.0f as a hard constraint). Both x and x^2 are exact in float,
so all the error is Horner rounding. Degrees 11 and 13 landed on the same
worst point (k = 61) with the same error, so degree 9 is where extra terms
stop paying. A Taylor series is the same evaluation with worse coefficients.

```c
int   k = ... folded, then negated as an integer when bit 7 is set (sin(pi) = +0)
float x = (float)k * (1.0f / 128.0f);    /* exact */
float t = x * x;                          /* exact: k*k <= 4096 */
float q =    0.0775605664f;               /* 0x1.3db026p-4 */
q = q * t - 0.598242283f;                 /* -0x1.324cd0p-1 */
q = q * t + 2.55006814f;                  /* 0x1.4668a2p+1 */
q = q * t - 5.16770983f;                  /* -0x1.4abbc2p+2 */
q = q * t + 3.14159274f;                  /* 0x1.921fb6p+1 */
return x * q;
```

Ways to get the polynomial to correct rounding that were not tried: evaluate
in double and round once (free in scalar x86, halves SIMD lane count), or a
fixed-point integer Horner in k with the coefficients tuned so every one of
the 65 results converts to the correctly rounded float.

### Two-level table

The angle-addition identity used as compression. The coarse table holds
sin(8 hi) for hi = 0..8; because cos(8 hi) = sin(64 - 8 hi) the same nine
entries serve cos read from the other end. The fine table holds sin and cos
of lo = 0..7. Multiples of 8 stay exact since the fine term vanishes there.
The quarter-wave fold is the last exact symmetry: anything under 65 entries
costs arithmetic and at least one rounding.

```c
static const float coarse[9] = {   /* sin(8 hi * pi/128) */
    0.0f, 0.195090324f, 0.382683426f, 0.555570245f, 0.707106769f,
    0.831469595f, 0.923879504f, 0.980785251f, 1.0f };
static const float fine[16] = {    /* sin(lo * pi/128) in 0..7, cos in 8..15 */
    0.0f, 0.024541229f, 0.0490676761f, 0.0735645667f,
    0.0980171412f, 0.122410677f, 0.146730468f, 0.170961887f,
    1.0f, 0.999698818f, 0.99879545f, 0.997290432f,
    0.99518472f, 0.992479563f, 0.989176512f, 0.985277653f };

unsigned hi = k >> 3, lo = k & 7u;
float s = coarse[hi] * fine[8u + lo] + coarse[8u - hi] * fine[lo];
```

All 25 values are entries of the existing 65-entry table (indices 0..8 by
eights, 0..7, and 64 down to 57).

## Sixteen at a time (AVX-512)

Same three formulations as 16-lane kernels, plus the 65-entry table through a
memory gather. Input was 16 raw bytes widened to 32-bit lanes with
`_mm512_cvtepu8_epi32`; the fold used `_mm512_test_epi32_mask` on bits 6 and
7 with `_mm512_mask_sub_epi32` for the mirror and `_mm512_mask_sub_ps` from
zero for the sign.

| kernel | ns/angle | constants in registers | worst error |
|---|---|---|---|
| scalar table, for reference | 0.69 | | 0.5 ulp |
| gather | 0.27 | 0 | 0.5 ulp |
| permute table | 0.10 | 4 zmm | 0.5 ulp, bit-identical to scalar |
| split | 0.14 | 2 zmm | 0.97 ulp |
| poly | 0.14 | 5 broadcasts | 1.29 ulp |

Ryzen 7 7700X, gcc 9, `-O2 -march=native`, 4096 angles resident in L1.

The register-resident table was the fastest and the most accurate. Entries
0..63 sit in four zmm registers; `_mm512_permutex2var_ps` indexes 32 entries
by bits 0..4 of k, so two of them and a blend on bit 5 cover 0..63, and entry
64 (= 1.0f) is a masked move from a broadcast. In a loop gcc hoists the four
table loads, so the table never touches memory again. The gather it replaces
costs more than the whole polynomial.

```c
__m512 t0 = _mm512_loadu_ps(table + 0),  t1 = _mm512_loadu_ps(table + 16);
__m512 t2 = _mm512_loadu_ps(table + 32), t3 = _mm512_loadu_ps(table + 48);
__m512 lo = _mm512_permutex2var_ps(t0, k, t1);
__m512 hi = _mm512_permutex2var_ps(t2, k, t3);
__m512 s  = _mm512_mask_blend_ps(_mm512_test_epi32_mask(k, set1(32)), lo, hi);
s = _mm512_mask_mov_ps(s, _mm512_cmpeq_epi32_mask(k, set1(64)), _mm512_set1_ps(1.0f));
```

The split kernel used one `_mm512_permutexvar_ps` each for sin(hi), cos(hi)
(index 8 - hi), sin(lo) and cos(lo) (index lo | 8), then `fmadd`. It tied the
polynomial on speed with better accuracy, so it would be the pick when four
registers are too many to give up.

Not tried: AVX2, where the permute would need eight ymm registers and the
gather may win instead.

## Testing that was in place

Exhaustive over all 256 angles for every scalar and vector variant: ulp error
against a reference that folds exactly (so sin(pi) compares against 0, not
against `sin(M_PI)` = 1.2e-16), bit-exact odd symmetry, half-turn negation,
even symmetry of cos, and the signed-zero cases at 0, pi, +-pi/2. Vector
gather and permute were checked bit-identical to the scalar table.
