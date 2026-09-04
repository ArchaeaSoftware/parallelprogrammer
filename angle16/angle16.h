/*
 * angle16.h -- header-only 16-bit angle type.
 *
 * The 16-bit edition of angle8.h. An angle is stored as a signed 16-bit
 * integer scaled by pi:
 *
 *     raw      radians              degrees
 *     -32768   -pi                  -180
 *     -16384   -pi/2                 -90
 *          0    0                      0
 *       8192    pi/4                  45
 *      16384    pi/2                  90
 *      32767    32767/32768 * pi     179.9945068359375
 *
 * 65536 raw units make one full turn, so two's-complement wrap around is
 * angular wrap around and angle addition, subtraction, negation and
 * "shortest signed difference" need no range checks. Resolution is
 * 1/65536 of a turn, 0.0055 degrees.
 *
 * Conversions to float are a single multiply. raw/32768 is a power-of-two
 * scale, so the constant pi/32768 has exactly the same mantissa as pi and
 * the scaling costs no precision. angle16_to_degrees() is exact: 180/32768
 * is 45 * 2^-13, and raw * 45 fits in a float mantissa.
 *
 * Bit tricks (see ANGLE16_TO_FLOAT_BITS / ANGLE16_FROM_FLOAT_BITS below):
 *
 *   int16 -> float: OR the sixteen bits into the mantissa of 2.0f so the
 *   float reads 2 + (raw+32768)/32768 = 3 + raw/32768, then subtract 3.0f.
 *   The subtraction is exact, so the result is bit-identical to
 *   (float)raw * (1/32768.0f) under every rounding mode.
 *
 *   float -> int16: add 1.5 * 2^23 so the FPU rounds the value to an integer
 *   in the low mantissa bits, then read the low sixteen bits. The mod-65536
 *   wrap comes for free. Rounds ties to even under the default rounding
 *   mode, which is also what llrintf() does, so the two paths agree. The
 *   trick covers |x| < 2^22, which is 128 turns; beyond that it falls back
 *   to the conversion.
 *
 * sin is a two-level table rather than the 16385-entry quarter wave (64 KB)
 * that the 8-bit edition's approach would need. After the quarter-wave fold
 * the index k = 128*hi + lo is split, and the angle-addition identity
 *
 *     sin(hi + lo) = sin(hi) + (cos(hi) sin(lo) - sin(hi) (1 - cos(lo)))
 *
 * is evaluated from a 129-entry coarse table of sin(hi * pi/256), which also
 * serves cos(hi) read from the other end, and two 128-entry fine tables of
 * sin(lo * pi/32768) and 1 - cos(lo * pi/32768). Storing 1 - cos keeps the
 * whole correction small, so the result is the coarse entry plus a small
 * term rounded once. 385 floats, 1540 bytes. The fold makes sin exactly
 * odd, cos exactly even, sin(pi) exactly +0 and sin(pi/2) exactly 1.
 *
 * Build options (define before including):
 *   ANGLE16_TO_FLOAT_BITS    1 (default) use the mantissa trick for int16->float
 *   ANGLE16_FROM_FLOAT_BITS  0 (default) use the magic-number trick for float->int16
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
#ifndef ANGLE16_H
#define ANGLE16_H

#include <stdint.h>
#include <string.h>
#include <math.h>

#ifndef ANGLE16_TO_FLOAT_BITS
#define ANGLE16_TO_FLOAT_BITS 1
#endif
#ifndef ANGLE16_FROM_FLOAT_BITS
#define ANGLE16_FROM_FLOAT_BITS 0
#endif

#ifndef ANGLE16_INLINE
#define ANGLE16_INLINE static inline
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------------ */
/* Type and constants                                                        */
/* ------------------------------------------------------------------------ */

typedef struct angle16_t { int16_t raw; } angle16_t;

#ifdef __cplusplus
#define ANGLE16_LIT(r) (angle16_t{ (int16_t)(r) })
#else
#define ANGLE16_LIT(r) ((angle16_t){ (int16_t)(r) })
#endif

#define ANGLE16_ZERO        ANGLE16_LIT(0)
#define ANGLE16_QUARTER_PI  ANGLE16_LIT(8192)
#define ANGLE16_HALF_PI     ANGLE16_LIT(16384)
#define ANGLE16_PI          ANGLE16_LIT(-32768)  /* +pi and -pi are the same angle */
#define ANGLE16_NEG_HALF_PI ANGLE16_LIT(-16384)

#define ANGLE16_PI_F        3.14159265358979323846f
#define ANGLE16_RAW_TO_RAD  (ANGLE16_PI_F / 32768.0f)  /* same mantissa as pi     */
#define ANGLE16_RAW_TO_DEG  0.0054931640625f           /* 180/32768, exact        */
#define ANGLE16_RAW_TO_TURN (1.0f / 65536.0f)
#define ANGLE16_RAW_TO_UNIT (1.0f / 32768.0f)          /* raw -> multiples of pi  */
#define ANGLE16_RAD_TO_UNIT (1.0f / ANGLE16_PI_F)
#define ANGLE16_DEG_TO_UNIT (1.0f / 180.0f)

/* ------------------------------------------------------------------------ */
/* Raw access                                                                */
/* ------------------------------------------------------------------------ */

ANGLE16_INLINE angle16_t angle16_from_raw(int16_t raw) { angle16_t a; a.raw = raw; return a; }
ANGLE16_INLINE int16_t   angle16_raw(angle16_t a)      { return a.raw; }
ANGLE16_INLINE uint16_t  angle16_uraw(angle16_t a)     { return (uint16_t)a.raw; } /* 0..65535, 0 = 0, 32768 = pi */

/* Reinterpret sixteen bits (0..65535 = 0..2pi) as an angle. */
ANGLE16_INLINE angle16_t angle16_from_uraw(uint16_t u)  { return angle16_from_raw((int16_t)u); }

/* ------------------------------------------------------------------------ */
/* int16 -> float                                                            */
/* ------------------------------------------------------------------------ */

/* raw/32768 in [-1, 1): the angle in multiples of pi. Plain conversion. */
ANGLE16_INLINE float angle16_to_unit_cast(angle16_t a)
{
    return (float)a.raw * ANGLE16_RAW_TO_UNIT;
}

/* raw/32768 in [-1, 1) via mantissa stuffing. Bit-identical to the cast form. */
ANGLE16_INLINE float angle16_to_unit_bits(angle16_t a)
{
    /* 0x40000000 is 2.0f. Mantissa bit 7 is worth 2 * 2^-16 = 1/32768 there,
     * so (raw ^ 0x8000) = raw + 32768 placed at bits 7..22 gives
     * 2 + (raw + 32768)/32768 = 3 + raw/32768. */
    uint32_t bits = 0x40000000u | ((uint32_t)((uint16_t)a.raw ^ 0x8000u) << 7);
    float f;
    memcpy(&f, &bits, sizeof f);
    return f - 3.0f;   /* exact: 3 + raw/32768 and raw/32768 are both representable */
}

ANGLE16_INLINE float angle16_to_unit(angle16_t a)
{
#if ANGLE16_TO_FLOAT_BITS
    return angle16_to_unit_bits(a);
#else
    return angle16_to_unit_cast(a);
#endif
}

ANGLE16_INLINE float angle16_to_radians(angle16_t a) { return angle16_to_unit(a) * ANGLE16_PI_F; }
ANGLE16_INLINE float angle16_to_degrees(angle16_t a) { return angle16_to_unit(a) * 180.0f; }   /* exact */
ANGLE16_INLINE float angle16_to_turns(angle16_t a)   { return angle16_to_unit(a) * 0.5f;   }   /* exact */

/* ------------------------------------------------------------------------ */
/* float -> int16                                                            */
/* ------------------------------------------------------------------------ */

/* Angle in multiples of pi -> raw, wrapped to [-32768, 32767], rounded to
 * nearest (ties to even). NaN and infinities map to 0. */
ANGLE16_INLINE angle16_t angle16_from_unit_cast(float u)
{
    float x = u * 32768.0f;
    /* Beyond 2^39 a float's ulp is >= 65536, so x is a multiple of a full
     * turn. The negated comparison also catches NaN. */
    if (!(fabsf(x) < 549755813888.0f)) return angle16_from_raw(0);
    return angle16_from_raw((int16_t)(uint16_t)llrintf(x));
}

/* Same result, via the magic-number rounding trick. Assumes the FPU is in
 * round-to-nearest mode (the default). */
ANGLE16_INLINE angle16_t angle16_from_unit_bits(float u)
{
    float x = u * 32768.0f;
    uint32_t bits;
    /* The trick needs |x| < 2^22 so that x + 1.5*2^23 stays in [2^23, 2^24)
     * where the float ulp is exactly 1. */
    if (!(fabsf(x) < 4194304.0f)) return angle16_from_unit_cast(u);
    x += 12582912.0f;                 /* 1.5 * 2^23: rounds x to an integer in the mantissa */
    memcpy(&bits, &x, sizeof bits);
    return angle16_from_raw((int16_t)(uint16_t)bits);   /* low 16 bits = round(x) mod 65536 */
}

ANGLE16_INLINE angle16_t angle16_from_unit(float u)
{
#if ANGLE16_FROM_FLOAT_BITS
    return angle16_from_unit_bits(u);
#else
    return angle16_from_unit_cast(u);
#endif
}

ANGLE16_INLINE angle16_t angle16_from_radians(float r) { return angle16_from_unit(r * ANGLE16_RAD_TO_UNIT); }
ANGLE16_INLINE angle16_t angle16_from_degrees(float d) { return angle16_from_unit(d * ANGLE16_DEG_TO_UNIT); }
ANGLE16_INLINE angle16_t angle16_from_turns(float t)   { return angle16_from_unit(t * 2.0f); }

/* Integer degrees -> raw without touching the FPU. Rounds half away from zero. */
ANGLE16_INLINE angle16_t angle16_from_degrees_int(int deg)
{
    int n;
    deg %= 360;                              /* keep deg * 8192 far from overflow */
    n = deg * 8192;                          /* raw = deg * 32768/180 = deg * 8192/45 */
    n = (n >= 0) ? (n + 22) / 45 : -((-n + 22) / 45);
    return angle16_from_raw((int16_t)(uint16_t)n);
}

/* ------------------------------------------------------------------------ */
/* Arithmetic (all mod 2pi, computed in unsigned to avoid signed overflow)   */
/* ------------------------------------------------------------------------ */

ANGLE16_INLINE angle16_t angle16_add(angle16_t a, angle16_t b)
{
    return angle16_from_raw((int16_t)(uint16_t)((uint16_t)a.raw + (uint16_t)b.raw));
}

/* a - b: the shortest signed rotation from b to a, in [-pi, pi). */
ANGLE16_INLINE angle16_t angle16_sub(angle16_t a, angle16_t b)
{
    return angle16_from_raw((int16_t)(uint16_t)((uint16_t)a.raw - (uint16_t)b.raw));
}

ANGLE16_INLINE angle16_t angle16_neg(angle16_t a)
{
    return angle16_from_raw((int16_t)(uint16_t)(0u - (uint16_t)a.raw));
}

/* Integer multiple of an angle, mod 2pi. */
ANGLE16_INLINE angle16_t angle16_scale(angle16_t a, int k)
{
    return angle16_from_raw((int16_t)(uint16_t)((unsigned)(uint16_t)a.raw * (unsigned)k));
}

ANGLE16_INLINE int angle16_eq(angle16_t a, angle16_t b) { return a.raw == b.raw; }

/* Unsigned distance from a to b going counter-clockwise, 0..65535. */
ANGLE16_INLINE uint16_t angle16_ccw_distance(angle16_t from, angle16_t to)
{
    return (uint16_t)((uint16_t)to.raw - (uint16_t)from.raw);
}

/* Absolute angular distance, 0..32768 (32768 = pi). */
ANGLE16_INLINE uint16_t angle16_distance(angle16_t a, angle16_t b)
{
    uint16_t d = angle16_ccw_distance(a, b);
    return (uint16_t)(d > 32768u ? 65536u - d : d);
}

/* Interpolate along the shortest arc. t is in 1/65536ths: 0 -> a, 65536 -> b.
 * The product reaches 2^31 in magnitude, so it is formed in 64 bits. */
ANGLE16_INLINE angle16_t angle16_lerp(angle16_t a, angle16_t b, int t)
{
    int d = angle16_sub(b, a).raw;          /* -32768..32767, shortest signed path */
    int64_t step = ((int64_t)d * t) >> 16;  /* arithmetic shift, rounds toward -inf */
    return angle16_add(a, angle16_from_raw((int16_t)(uint16_t)(uint64_t)step));
}

/* ------------------------------------------------------------------------ */
/* Trigonometry via a two-level quarter-wave table                           */
/* ------------------------------------------------------------------------ */

/* Quarter-wave reduction: sin(u * pi/32768) = +-sin(k * pi/32768) with k in
 * 0..16384. Bit 14 of u mirrors the index, bit 15 flips the sign. */
ANGLE16_INLINE unsigned angle16__quarter(uint16_t u)
{
    unsigned i = u & 16383u;
    return (u & 16384u) ? 16384u - i : i;
}

/* sin(hi * pi/256) for hi = 0..128, rounded to float. cos(hi * pi/256) is
 * the same table read from the other end. */
static const float angle16__sin_coarse[129] = {
    0.0f, 0.0122715384f, 0.024541229f, 0.0368072242f,
    0.0490676761f, 0.061320737f, 0.0735645667f, 0.0857973099f,
    0.0980171412f, 0.110222206f, 0.122410677f, 0.134580702f,
    0.146730468f, 0.15885815f, 0.170961887f, 0.183039889f,
    0.195090324f, 0.207111374f, 0.219101235f, 0.231058106f,
    0.242980182f, 0.254865646f, 0.266712755f, 0.27851969f,
    0.290284663f, 0.302005947f, 0.313681751f, 0.32531029f,
    0.336889863f, 0.348418683f, 0.359895051f, 0.371317208f,
    0.382683426f, 0.393992037f, 0.405241311f, 0.416429549f,
    0.427555084f, 0.438616246f, 0.449611336f, 0.460538715f,
    0.471396744f, 0.482183784f, 0.492898196f, 0.50353837f,
    0.514102757f, 0.524589658f, 0.534997642f, 0.545324981f,
    0.555570245f, 0.565731823f, 0.575808167f, 0.585797846f,
    0.59569931f, 0.605511069f, 0.615231574f, 0.624859512f,
    0.634393275f, 0.643831551f, 0.653172851f, 0.662415802f,
    0.671558976f, 0.680601001f, 0.689540565f, 0.698376238f,
    0.707106769f, 0.715730846f, 0.724247098f, 0.732654274f,
    0.740951121f, 0.749136388f, 0.757208824f, 0.765167236f,
    0.773010433f, 0.780737221f, 0.78834641f, 0.795836926f,
    0.803207517f, 0.81045717f, 0.817584813f, 0.824589312f,
    0.831469595f, 0.838224709f, 0.84485358f, 0.851355195f,
    0.857728601f, 0.863972843f, 0.870086968f, 0.876070082f,
    0.881921291f, 0.887639642f, 0.893224299f, 0.898674488f,
    0.903989315f, 0.909168005f, 0.914209783f, 0.919113874f,
    0.923879504f, 0.928506076f, 0.932992816f, 0.937339008f,
    0.941544056f, 0.945607305f, 0.949528158f, 0.953306019f,
    0.956940353f, 0.960430503f, 0.963776052f, 0.966976464f,
    0.970031261f, 0.972939968f, 0.975702107f, 0.97831738f,
    0.980785251f, 0.983105481f, 0.985277653f, 0.987301409f,
    0.989176512f, 0.990902662f, 0.992479563f, 0.993906975f,
    0.99518472f, 0.996312618f, 0.997290432f, 0.998118103f,
    0.99879545f, 0.999322355f, 0.999698818f, 0.999924719f,
    1.0f,
};

/* sin(lo * pi/32768) for lo = 0..127, rounded to float. */
static const float angle16__sin_fine[128] = {
    0.0f, 9.58738019e-05f, 0.000191747604f, 0.000287621398f,
    0.000383495179f, 0.000479368988f, 0.000575242739f, 0.000671116519f,
    0.000766990299f, 0.000862864079f, 0.000958737859f, 0.00105461164f,
    0.00115048536f, 0.00124635908f, 0.0013422328f, 0.00143810653f,
    0.00153398013f, 0.00162985385f, 0.00172572758f, 0.00182160118f,
    0.00191747479f, 0.00201334851f, 0.00210922211f, 0.00220509549f,
    0.00230096909f, 0.0023968427f, 0.0024927163f, 0.00258858968f,
    0.00268446305f, 0.00278033665f, 0.00287621003f, 0.0029720834f,
    0.00306795677f, 0.00316383014f, 0.00325970352f, 0.00335557666f,
    0.00345145003f, 0.00354732317f, 0.00364319631f, 0.00373906945f,
    0.00383494259f, 0.0039308155f, 0.00402668864f, 0.00412256178f,
    0.00421843445f, 0.00431430759f, 0.00441018026f, 0.0045060534f,
    0.00460192608f, 0.00469779875f, 0.00479367143f, 0.0048895441f,
    0.00498541677f, 0.00508128945f, 0.00517716212f, 0.00527303433f,
    0.005368907f, 0.00546477921f, 0.00556065189f, 0.0056565241f,
    0.0057523963f, 0.00584826851f, 0.00594414072f, 0.00604001246f,
    0.00613588467f, 0.00623175642f, 0.00632762862f, 0.00642350037f,
    0.00651937211f, 0.00661524385f, 0.00671111559f, 0.00680698734f,
    0.00690285861f, 0.00699873036f, 0.00709460163f, 0.00719047291f,
    0.00728634419f, 0.00738221547f, 0.00747808674f, 0.00757395755f,
    0.00766982883f, 0.00776569964f, 0.00786157046f, 0.00795744173f,
    0.00805331208f, 0.00814918242f, 0.0082450537f, 0.00834092405f,
    0.00843679439f, 0.00853266474f, 0.00862853508f, 0.00872440543f,
    0.00882027484f, 0.00891614519f, 0.00901201554f, 0.00910788495f,
    0.00920375437f, 0.00929962471f, 0.00939549413f, 0.00949136354f,
    0.00958723295f, 0.00968310237f, 0.00977897178f, 0.0098748412f,
    0.00997070968f, 0.0100665791f, 0.0101624476f, 0.010258317f,
    0.0103541855f, 0.010450054f, 0.0105459224f, 0.0106417909f,
    0.0107376594f, 0.010833527f, 0.0109293954f, 0.0110252639f,
    0.0111211315f, 0.011216999f, 0.0113128666f, 0.0114087351f,
    0.0115046017f, 0.0116004692f, 0.0116963368f, 0.0117922043f,
    0.011888071f, 0.0119839376f, 0.0120798051f, 0.0121756718f,
};

/* 1 - cos(lo * pi/32768) for lo = 0..127, rounded to float. At most 7.5e-5,
 * so the rounding of the stored value is far below a float ulp of any sine
 * it is multiplied into. */
static const float angle16__cm1_fine[128] = {
    0.0f, 4.59589256e-09f, 1.83835702e-08f, 4.13630339e-08f,
    7.35342809e-08f, 1.14897318e-07f, 1.65452136e-07f, 2.25198733e-07f,
    2.94137124e-07f, 3.72267294e-07f, 4.59589245e-07f, 5.56102975e-07f,
    6.61808485e-07f, 7.76705747e-07f, 9.00794817e-07f, 1.03407569e-06f,
    1.17654827e-06f, 1.32821265e-06f, 1.48906884e-06f, 1.65911683e-06f,
    1.83835652e-06f, 2.02678802e-06f, 2.22441122e-06f, 2.43122622e-06f,
    2.64723303e-06f, 2.87243165e-06f, 3.10682185e-06f, 3.35040386e-06f,
    3.60317767e-06f, 3.8651433e-06f, 4.1363005e-06f, 4.41664952e-06f,
    4.70619034e-06f, 5.00492297e-06f, 5.31284741e-06f, 5.6299632e-06f,
    5.9562708e-06f, 6.29177066e-06f, 6.63646188e-06f, 6.99034445e-06f,
    7.35341928e-06f, 7.72568546e-06f, 8.10714391e-06f, 8.49779371e-06f,
    8.89763487e-06f, 9.30666829e-06f, 9.72489306e-06f, 1.01523101e-05f,
    1.05889185e-05f, 1.10347182e-05f, 1.14897093e-05f, 1.19538927e-05f,
    1.24272683e-05f, 1.29098344e-05f, 1.34015936e-05f, 1.39025433e-05f,
    1.44126852e-05f, 1.49320185e-05f, 1.54605423e-05f, 1.59982592e-05f,
    1.65451675e-05f, 1.71012671e-05f, 1.76665599e-05f, 1.82410422e-05f,
    1.88247177e-05f, 1.94175846e-05f, 2.00196409e-05f, 2.06308905e-05f,
    2.12513332e-05f, 2.18809655e-05f, 2.25197891e-05f, 2.31678059e-05f,
    2.38250122e-05f, 2.44914118e-05f, 2.51670026e-05f, 2.58517848e-05f,
    2.65457584e-05f, 2.72489233e-05f, 2.79612814e-05f, 2.8682829e-05f,
    2.94135698e-05f, 3.01535001e-05f, 3.09026218e-05f, 3.16609367e-05f,
    3.24284447e-05f, 3.32051422e-05f, 3.39910293e-05f, 3.47861096e-05f,
    3.5590383e-05f, 3.6403846e-05f, 3.72264985e-05f, 3.80583442e-05f,
    3.88993831e-05f, 3.97496115e-05f, 4.0609033e-05f, 4.14776441e-05f,
    4.23554484e-05f, 4.32424422e-05f, 4.41386292e-05f, 4.50440057e-05f,
    4.59585754e-05f, 4.68823346e-05f, 4.7815287e-05f, 4.87574289e-05f,
    4.9708764e-05f, 5.06692886e-05f, 5.16390064e-05f, 5.26179138e-05f,
    5.36060143e-05f, 5.46033043e-05f, 5.56097875e-05f, 5.66254603e-05f,
    5.76503226e-05f, 5.8684378e-05f, 5.97276266e-05f, 6.07800648e-05f,
    6.18416962e-05f, 6.29125134e-05f, 6.39925274e-05f, 6.5081731e-05f,
    6.61801241e-05f, 6.72877068e-05f, 6.84044862e-05f, 6.95304552e-05f,
    7.06656137e-05f, 7.18099618e-05f, 7.29635067e-05f, 7.41262338e-05f,
};

ANGLE16_INLINE float angle16_sin(angle16_t a)
{
    uint16_t u = (uint16_t)a.raw;           /* 0..65535 = 0..2pi */
    unsigned k = angle16__quarter(u), hi = k >> 7, lo = k & 127u;
    float s_hi = angle16__sin_coarse[hi];
    float c_hi = angle16__sin_coarse[128u - hi];
    float s = s_hi + (c_hi * angle16__sin_fine[lo] - s_hi * angle16__cm1_fine[lo]);
    return (u & 32768u) ? 0.0f - s : s;     /* 0.0f - s avoids producing -0.0 */
}

ANGLE16_INLINE float angle16_cos(angle16_t a)
{
    return angle16_sin(angle16_from_uraw((uint16_t)((uint16_t)a.raw + 16384u)));
}

/* Nearest 16-bit angle to atan2(y, x). Uses libm's atan2f. */
ANGLE16_INLINE angle16_t angle16_atan2(float y, float x)
{
    return angle16_from_radians(atan2f(y, x));
}

#ifdef __cplusplus
} /* extern "C" */

inline angle16_t  operator+ (angle16_t a, angle16_t b) { return angle16_add(a, b); }
inline angle16_t  operator- (angle16_t a, angle16_t b) { return angle16_sub(a, b); }
inline angle16_t  operator- (angle16_t a)              { return angle16_neg(a); }
inline angle16_t  operator* (angle16_t a, int k)       { return angle16_scale(a, k); }
inline angle16_t  operator* (int k, angle16_t a)       { return angle16_scale(a, k); }
inline angle16_t& operator+=(angle16_t& a, angle16_t b) { a = angle16_add(a, b); return a; }
inline angle16_t& operator-=(angle16_t& a, angle16_t b) { a = angle16_sub(a, b); return a; }
inline angle16_t& operator*=(angle16_t& a, int k)       { a = angle16_scale(a, k); return a; }
inline bool       operator==(angle16_t a, angle16_t b) { return a.raw == b.raw; }
inline bool       operator!=(angle16_t a, angle16_t b) { return a.raw != b.raw; }
#endif

#endif /* ANGLE16_H */
