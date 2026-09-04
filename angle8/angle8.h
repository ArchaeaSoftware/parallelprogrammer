/*
 * angle8.h -- header-only 8-bit angle type.
 *
 * An angle is stored as a signed 8-bit integer scaled by pi:
 *
 *     raw    radians          degrees
 *     -128   -pi              -180
 *      -64   -pi/2             -90
 *        0    0                  0
 *       32    pi/4              45
 *       64    pi/2              90
 *      127    127/128 * pi     178.59375
 *
 * Because 256 raw units make one full turn, ordinary two's-complement wrap
 * around IS angular wrap around: adding 1 to 127 (just below +pi) yields
 * -128 (-pi), which is the same direction. Angle addition, subtraction,
 * negation and "shortest signed difference" therefore need no range checks.
 *
 * Conversions to float are a single multiply. raw/128 is a power-of-two
 * scale, so the constant pi/128 has exactly the same mantissa as pi and
 * the scaling costs no precision. angle8_to_degrees() is exact: 180/128 =
 * 1.40625 is representable.
 *
 * Bit tricks (see ANGLE8_TO_FLOAT_BITS / ANGLE8_FROM_FLOAT_BITS below):
 *
 *   int8 -> float: OR the byte into the mantissa of 2.0f so the float reads
 *   2 + (raw+128)/128 = 3 + raw/128, then subtract 3.0f. The subtraction is
 *   exact, so the result is bit-identical to (float)raw * (1/128.0f) under
 *   every rounding mode. Useful on targets where int->float conversion is
 *   slow or absent.
 *
 *   float -> int8: add 1.5 * 2^23 so the FPU rounds the value to an integer
 *   in the low mantissa bits, then read the low byte. The mod-256 wrap comes
 *   for free from taking the low byte. Rounds ties to even under the default
 *   rounding mode, which is also what lrintf() does, so the two paths agree.
 *
 * Build options (define before including):
 *   ANGLE8_TO_FLOAT_BITS    1 (default) use the mantissa trick for int8->float
 *   ANGLE8_FROM_FLOAT_BITS  0 (default) use the magic-number trick for float->int8
 *
 * C99 or C++11. Requires <math.h> (link -lm on some C toolchains).
 *
 * Copyright (C) 2026 by Archaea Software, LLC.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice, this
 *    list of conditions and the following disclaimer in the documentation and/or
 *    other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef ANGLE8_H
#define ANGLE8_H

#include <stdint.h>
#include <string.h>
#include <math.h>

#ifndef ANGLE8_TO_FLOAT_BITS
#define ANGLE8_TO_FLOAT_BITS 1
#endif
#ifndef ANGLE8_FROM_FLOAT_BITS
#define ANGLE8_FROM_FLOAT_BITS 0
#endif

#ifndef ANGLE8_INLINE
#define ANGLE8_INLINE static inline
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------------ */
/* Type and constants                                                        */
/* ------------------------------------------------------------------------ */

typedef struct angle8_t { int8_t raw; } angle8_t;

#ifdef __cplusplus
#define ANGLE8_LIT(r) (angle8_t{ (int8_t)(r) })
#else
#define ANGLE8_LIT(r) ((angle8_t){ (int8_t)(r) })
#endif

#define ANGLE8_ZERO        ANGLE8_LIT(0)
#define ANGLE8_QUARTER_PI  ANGLE8_LIT(32)
#define ANGLE8_HALF_PI     ANGLE8_LIT(64)
#define ANGLE8_PI          ANGLE8_LIT(-128)   /* +pi and -pi are the same angle */
#define ANGLE8_NEG_HALF_PI ANGLE8_LIT(-64)

#define ANGLE8_PI_F        3.14159265358979323846f
#define ANGLE8_RAW_TO_RAD  (ANGLE8_PI_F / 128.0f)   /* same mantissa as pi     */
#define ANGLE8_RAW_TO_DEG  1.40625f                 /* 180/128, exact          */
#define ANGLE8_RAW_TO_TURN (1.0f / 256.0f)
#define ANGLE8_RAW_TO_UNIT (1.0f / 128.0f)          /* raw -> multiples of pi  */
#define ANGLE8_RAD_TO_UNIT (1.0f / ANGLE8_PI_F)
#define ANGLE8_DEG_TO_UNIT (1.0f / 180.0f)

/* ------------------------------------------------------------------------ */
/* Raw access                                                                */
/* ------------------------------------------------------------------------ */

ANGLE8_INLINE angle8_t angle8_from_raw(int8_t raw) { angle8_t a; a.raw = raw; return a; }
ANGLE8_INLINE int8_t   angle8_raw(angle8_t a)      { return a.raw; }
ANGLE8_INLINE uint8_t  angle8_uraw(angle8_t a)     { return (uint8_t)a.raw; } /* 0..255, 0 = 0, 128 = pi */

/* Reinterpret a byte (0..255 = 0..2pi) as an angle. */
ANGLE8_INLINE angle8_t angle8_from_uraw(uint8_t u)  { return angle8_from_raw((int8_t)u); }

/* ------------------------------------------------------------------------ */
/* int8 -> float                                                             */
/* ------------------------------------------------------------------------ */

/* raw/128 in [-1, 1): the angle in multiples of pi. Plain conversion. */
ANGLE8_INLINE float angle8_to_unit_cast(angle8_t a)
{
    return (float)a.raw * ANGLE8_RAW_TO_UNIT;
}

/* raw/128 in [-1, 1) via mantissa stuffing. Bit-identical to the cast form. */
ANGLE8_INLINE float angle8_to_unit_bits(angle8_t a)
{
    /* 0x40000000 is 2.0f. Mantissa bit 15 is worth 2 * 2^-8 = 1/128 there,
     * so byte (raw ^ 0x80) = raw + 128 placed at bits 15..22 gives
     * 2 + (raw + 128)/128 = 3 + raw/128. */
    uint32_t bits = 0x40000000u | ((uint32_t)((uint8_t)a.raw ^ 0x80u) << 15);
    float f;
    memcpy(&f, &bits, sizeof f);
    return f - 3.0f;   /* exact: 3 + raw/128 and raw/128 are both representable */
}

ANGLE8_INLINE float angle8_to_unit(angle8_t a)
{
#if ANGLE8_TO_FLOAT_BITS
    return angle8_to_unit_bits(a);
#else
    return angle8_to_unit_cast(a);
#endif
}

ANGLE8_INLINE float angle8_to_radians(angle8_t a) { return angle8_to_unit(a) * ANGLE8_PI_F; }
ANGLE8_INLINE float angle8_to_degrees(angle8_t a) { return angle8_to_unit(a) * 180.0f; }   /* exact */
ANGLE8_INLINE float angle8_to_turns(angle8_t a)   { return angle8_to_unit(a) * 0.5f;   }   /* exact */

/* ------------------------------------------------------------------------ */
/* float -> int8                                                             */
/* ------------------------------------------------------------------------ */

/* Angle in multiples of pi -> raw, wrapped to [-128, 127], rounded to nearest
 * (ties to even). NaN and infinities map to 0. */
ANGLE8_INLINE angle8_t angle8_from_unit_cast(float u)
{
    float x = u * 128.0f;
    /* Beyond 2^31 a float's ulp is >= 256, so x is a multiple of a full turn.
     * The negated comparison also catches NaN. */
    if (!(fabsf(x) < 2147483648.0f)) return angle8_from_raw(0);
    return angle8_from_raw((int8_t)(uint8_t)lrintf(x));
}

/* Same result, via the magic-number rounding trick. Assumes the FPU is in
 * round-to-nearest mode (the default). */
ANGLE8_INLINE angle8_t angle8_from_unit_bits(float u)
{
    float x = u * 128.0f;
    uint32_t bits;
    /* The trick needs |x| < 2^22 so that x + 1.5*2^23 stays in [2^23, 2^24)
     * where the float ulp is exactly 1. */
    if (!(fabsf(x) < 4194304.0f)) return angle8_from_unit_cast(u);
    x += 12582912.0f;                 /* 1.5 * 2^23: rounds x to an integer in the mantissa */
    memcpy(&bits, &x, sizeof bits);
    return angle8_from_raw((int8_t)(uint8_t)bits);   /* low byte = round(x) mod 256 */
}

ANGLE8_INLINE angle8_t angle8_from_unit(float u)
{
#if ANGLE8_FROM_FLOAT_BITS
    return angle8_from_unit_bits(u);
#else
    return angle8_from_unit_cast(u);
#endif
}

ANGLE8_INLINE angle8_t angle8_from_radians(float r) { return angle8_from_unit(r * ANGLE8_RAD_TO_UNIT); }
ANGLE8_INLINE angle8_t angle8_from_degrees(float d) { return angle8_from_unit(d * ANGLE8_DEG_TO_UNIT); }
ANGLE8_INLINE angle8_t angle8_from_turns(float t)   { return angle8_from_unit(t * 2.0f); }

/* Integer degrees -> raw without touching the FPU. Rounds half away from zero. */
ANGLE8_INLINE angle8_t angle8_from_degrees_int(int deg)
{
    int n;
    deg %= 360;                              /* keep deg * 32 far from overflow */
    n = deg * 32;                            /* raw = deg * 128/180 = deg * 32/45 */
    n = (n >= 0) ? (n + 22) / 45 : -((-n + 22) / 45);
    return angle8_from_raw((int8_t)(uint8_t)n);
}

/* ------------------------------------------------------------------------ */
/* Arithmetic (all mod 2pi, computed in unsigned to avoid signed overflow)   */
/* ------------------------------------------------------------------------ */

ANGLE8_INLINE angle8_t angle8_add(angle8_t a, angle8_t b)
{
    return angle8_from_raw((int8_t)(uint8_t)((uint8_t)a.raw + (uint8_t)b.raw));
}

/* a - b: the shortest signed rotation from b to a, in [-pi, pi). */
ANGLE8_INLINE angle8_t angle8_sub(angle8_t a, angle8_t b)
{
    return angle8_from_raw((int8_t)(uint8_t)((uint8_t)a.raw - (uint8_t)b.raw));
}

ANGLE8_INLINE angle8_t angle8_neg(angle8_t a)
{
    return angle8_from_raw((int8_t)(uint8_t)(0u - (uint8_t)a.raw));
}

/* Integer multiple of an angle, mod 2pi. */
ANGLE8_INLINE angle8_t angle8_scale(angle8_t a, int k)
{
    return angle8_from_raw((int8_t)(uint8_t)((unsigned)(uint8_t)a.raw * (unsigned)k));
}

ANGLE8_INLINE int angle8_eq(angle8_t a, angle8_t b) { return a.raw == b.raw; }

/* Unsigned distance from a to b going counter-clockwise, 0..255. */
ANGLE8_INLINE uint8_t angle8_ccw_distance(angle8_t from, angle8_t to)
{
    return (uint8_t)((uint8_t)to.raw - (uint8_t)from.raw);
}

/* Absolute angular distance, 0..128 (128 = pi). */
ANGLE8_INLINE uint8_t angle8_distance(angle8_t a, angle8_t b)
{
    uint8_t d = angle8_ccw_distance(a, b);
    return (uint8_t)(d > 128u ? 256u - d : d);
}

/* Interpolate along the shortest arc. t is in 1/256ths: 0 -> a, 256 -> b. */
ANGLE8_INLINE angle8_t angle8_lerp(angle8_t a, angle8_t b, int t)
{
    int d = angle8_sub(b, a).raw;           /* -128..127, shortest signed path */
    int step = (d * t) >> 8;                /* arithmetic shift, rounds toward -inf */
    return angle8_add(a, angle8_from_raw((int8_t)(uint8_t)(unsigned)step));
}

/* ------------------------------------------------------------------------ */
/* Trigonometry via a 65-entry quarter-wave table                            */
/* ------------------------------------------------------------------------ */

/* sin(k * pi/128) for k = 0..64, rounded to float. */
static const float angle8__sin_quarter[65] = {
    0.0f, 0.024541229009628296f, 0.049067676067352295f, 0.0735645666718483f,
    0.0980171412229538f, 0.12241067737340927f, 0.1467304676771164f, 0.1709618866443634f,
    0.19509032368659973f, 0.21910123527050018f, 0.24298018217086792f, 0.2667127549648285f,
    0.290284663438797f, 0.3136817514896393f, 0.3368898630142212f, 0.3598950505256653f,
    0.3826834261417389f, 0.40524131059646606f, 0.4275550842285156f, 0.4496113359928131f,
    0.4713967442512512f, 0.49289819598197937f, 0.5141027569770813f, 0.5349976420402527f,
    0.5555702447891235f, 0.5758081674575806f, 0.5956993103027344f, 0.6152315735816956f,
    0.6343932747840881f, 0.6531728506088257f, 0.6715589761734009f, 0.6895405650138855f,
    0.7071067690849304f, 0.7242470979690552f, 0.7409511208534241f, 0.7572088241577148f,
    0.7730104327201843f, 0.7883464097976685f, 0.803207516670227f, 0.8175848126411438f,
    0.8314695954322815f, 0.8448535799980164f, 0.8577286005020142f, 0.8700869679450989f,
    0.8819212913513184f, 0.89322429895401f, 0.903989315032959f, 0.91420978307724f,
    0.9238795042037964f, 0.9329928159713745f, 0.9415440559387207f, 0.949528157711029f,
    0.9569403529167175f, 0.9637760519981384f, 0.9700312614440918f, 0.9757021069526672f,
    0.9807852506637573f, 0.9852776527404785f, 0.9891765117645264f, 0.9924795627593994f,
    0.9951847195625305f, 0.9972904324531555f, 0.9987954497337341f, 0.99969881772995f,
    1.0f,
};

ANGLE8_INLINE float angle8_sin(angle8_t a)
{
    uint8_t u = (uint8_t)a.raw;             /* 0..255 = 0..2pi */
    uint8_t i = u & 63u;
    float s = angle8__sin_quarter[(u & 64u) ? (64u - i) : i];
    return (u & 128u) ? 0.0f - s : s;   /* 0.0f - s avoids producing -0.0 */
}

ANGLE8_INLINE float angle8_cos(angle8_t a)
{
    return angle8_sin(angle8_from_uraw((uint8_t)((uint8_t)a.raw + 64u)));
}

/* Nearest 8-bit angle to atan2(y, x). Uses libm's atan2f. */
ANGLE8_INLINE angle8_t angle8_atan2(float y, float x)
{
    return angle8_from_radians(atan2f(y, x));
}

#ifdef __cplusplus
} /* extern "C" */

inline angle8_t  operator+ (angle8_t a, angle8_t b) { return angle8_add(a, b); }
inline angle8_t  operator- (angle8_t a, angle8_t b) { return angle8_sub(a, b); }
inline angle8_t  operator- (angle8_t a)             { return angle8_neg(a); }
inline angle8_t  operator* (angle8_t a, int k)      { return angle8_scale(a, k); }
inline angle8_t  operator* (int k, angle8_t a)      { return angle8_scale(a, k); }
inline angle8_t& operator+=(angle8_t& a, angle8_t b) { a = angle8_add(a, b); return a; }
inline angle8_t& operator-=(angle8_t& a, angle8_t b) { a = angle8_sub(a, b); return a; }
inline angle8_t& operator*=(angle8_t& a, int k)      { a = angle8_scale(a, k); return a; }
inline bool      operator==(angle8_t a, angle8_t b) { return a.raw == b.raw; }
inline bool      operator!=(angle8_t a, angle8_t b) { return a.raw != b.raw; }
#endif

#endif /* ANGLE8_H */
