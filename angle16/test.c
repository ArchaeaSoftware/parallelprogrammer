/*
 * test.c
 *
 * Exhaustive tests for angle16.h. Build and run with `make test`.
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
#include "angle16.h"
#include <stdio.h>
#include <stdlib.h>

#define PI_D 3.14159265358979323846

static int fails = 0, checks = 0;
#define CHECK(cond) do { checks++; if (!(cond)) { fails++; \
    printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); } } while (0)

static int same_bits(float a, float b)
{
    uint32_t x, y;
    if (a == 0.0f && b == 0.0f) return 1;   /* treat +0 and -0 alike */
    memcpy(&x, &a, 4); memcpy(&y, &b, 4);
    return x == y;
}
static angle16_t A(int r) { return angle16_from_raw((int16_t)r); }

static void test_to_float(void)
{
    int r;
    for (r = -32768; r <= 32767; r++) {
        angle16_t a = A(r);
        float unit = angle16_to_unit(a);
        double ref = r * PI_D / 32768.0;
        float rad = angle16_to_radians(a);
        CHECK(same_bits(angle16_to_unit_bits(a), angle16_to_unit_cast(a)));
        CHECK(unit == (float)r / 32768.0f);
        CHECK(fabs(rad - ref) <= fabs(ref) * 1.2e-7);   /* within 1 float ulp */
        CHECK(angle16_to_degrees(a) == (float)r * 0.0054931640625f);
        CHECK(angle16_to_turns(a) == (float)r / 65536.0f);
    }
    /* Exact endpoints and landmarks. */
    CHECK(angle16_to_radians(A(-32768)) == -ANGLE16_PI_F);
    CHECK(angle16_to_radians(A(32767)) == (32767.0f / 32768.0f) * ANGLE16_PI_F);
    CHECK(angle16_to_radians(A(16384)) == ANGLE16_PI_F * 0.5f);
    CHECK(angle16_to_radians(A(0)) == 0.0f);
    CHECK(angle16_to_degrees(A(16384)) == 90.0f);
    CHECK(angle16_to_degrees(A(-16384)) == -90.0f);
    CHECK(angle16_to_degrees(A(-32768)) == -180.0f);
    CHECK(angle16_to_degrees(A(32767)) == 179.9945068359375f);
    CHECK(angle16_to_degrees(A(1)) == 0.0054931640625f);
    CHECK(angle16_to_turns(A(16384)) == 0.25f);
    CHECK(angle16_to_turns(A(-32768)) == -0.5f);
}

static void test_from_float(void)
{
    int r, j, k;
    uint32_t seed = 12345;

    /* Round trips through every representation, both implementations. */
    for (r = -32768; r <= 32767; r++) {
        angle16_t a = A(r);
        CHECK(angle16_from_radians(angle16_to_radians(a)).raw == r);
        CHECK(angle16_from_degrees(angle16_to_degrees(a)).raw == r);
        CHECK(angle16_from_turns(angle16_to_turns(a)).raw == r);
        CHECK(angle16_from_unit_cast(angle16_to_unit(a)).raw == r);
        CHECK(angle16_from_unit_bits(angle16_to_unit(a)).raw == r);
    }

    /* Bit-trick and cast paths agree exactly on exact ties (k/65536 -> x = k/2). */
    for (k = -200000; k <= 200000; k++) {
        float u = (float)k / 65536.0f;
        CHECK(angle16_from_unit_bits(u).raw == angle16_from_unit_cast(u).raw);
    }
    /* ...and on pseudo-random values across the trick's range and past it. */
    for (j = 0; j < 2000000; j++) {
        float u;
        seed = seed * 1664525u + 1013904223u;
        memcpy(&u, &seed, 4);
        if (!(fabsf(u) < 150.0f)) u = (float)(seed >> 8) * (1.0f / 16777216.0f) * 64.0f - 32.0f;
        CHECK(angle16_from_unit_bits(u).raw == angle16_from_unit_cast(u).raw);
    }

    /* Ties round to even. */
    CHECK(angle16_from_unit( 0.5f / 32768.0f).raw == 0);
    CHECK(angle16_from_unit( 1.5f / 32768.0f).raw == 2);
    CHECK(angle16_from_unit(-0.5f / 32768.0f).raw == 0);
    CHECK(angle16_from_unit(-1.5f / 32768.0f).raw == -2);

    /* Wrap: adding whole turns changes nothing (sums chosen to be exact). */
    for (j = -512; j <= 512; j++) {
        float t = (float)j / 1024.0f;
        angle16_t base = angle16_from_turns(t);
        for (k = -3; k <= 3; k++) {
            CHECK(angle16_from_turns(t + (float)k).raw == base.raw);
            CHECK(angle16_from_unit_bits((t + (float)k) * 2.0f).raw == base.raw);
        }
    }
    CHECK(angle16_from_unit(1.0f).raw == -32768);
    CHECK(angle16_from_unit(-1.0f).raw == -32768);
    CHECK(angle16_from_unit(2.0f).raw == 0);
    CHECK(angle16_from_unit(3.0f).raw == -32768);
    CHECK(angle16_from_degrees(360.0f).raw == 0);
    CHECK(angle16_from_degrees(270.0f).raw == -16384);
    CHECK(angle16_from_degrees(-90.0f).raw == -16384);
    CHECK(angle16_from_radians(ANGLE16_PI_F).raw == -32768);
    CHECK(angle16_from_radians(-ANGLE16_PI_F).raw == -32768);

    /* Non-finite and huge inputs are well defined. */
    CHECK(angle16_from_radians(NAN).raw == 0);
    CHECK(angle16_from_radians(INFINITY).raw == 0);
    CHECK(angle16_from_radians(-INFINITY).raw == 0);
    CHECK(angle16_from_unit_bits(NAN).raw == 0);
    CHECK(angle16_from_unit(1e30f).raw == 0);
    CHECK(angle16_from_unit(1e8f).raw == 0);                     /* |x| >= 2^39: multiple of 65536 */
    CHECK(angle16_from_unit(1e6f).raw == 0);                     /* x = 15625 * 2^21 */
    CHECK(angle16_from_unit(256.000030517578125f).raw == 1);     /* x = 2^23 + 1 */
    CHECK(angle16_from_unit_bits(256.000030517578125f).raw == 1);/* falls back past 2^22 */
    CHECK(angle16_from_unit_bits(1e7f).raw == angle16_from_unit_cast(1e7f).raw);
}

static void test_arithmetic(void)
{
    int i, j;
    CHECK(angle16_add(A(32767), A(1)).raw == -32768);
    CHECK(angle16_add(A(-32768), A(-1)).raw == 32767);
    CHECK(angle16_sub(A(-32768), A(1)).raw == 32767);
    CHECK(angle16_sub(A(32767), A(-1)).raw == -32768);
    CHECK(angle16_neg(A(-32768)).raw == -32768);
    CHECK(angle16_neg(A(0)).raw == 0);
    CHECK(angle16_neg(A(1)).raw == -1);
    CHECK(angle16_scale(A(16384), 2).raw == -32768);
    CHECK(angle16_scale(A(16384), 4).raw == 0);
    CHECK(angle16_scale(A(1), 65536).raw == 0);
    CHECK(angle16_scale(A(3), -1).raw == -3);
    CHECK(angle16_scale(A(-32768), -1).raw == -32768);
    CHECK(angle16_eq(ANGLE16_PI, A(-32768)));
    CHECK(ANGLE16_HALF_PI.raw == 16384 && ANGLE16_NEG_HALF_PI.raw == -16384);
    CHECK(ANGLE16_QUARTER_PI.raw == 8192 && ANGLE16_ZERO.raw == 0);

    /* Shortest signed difference crosses the wrap. */
    CHECK(angle16_sub(A(-30000), A(30000)).raw == 5536);
    CHECK(angle16_sub(A(30000), A(-30000)).raw == -5536);
    CHECK(angle16_ccw_distance(A(30000), A(-30000)) == 5536);
    CHECK(angle16_ccw_distance(A(-30000), A(30000)) == 60000);
    CHECK(angle16_distance(A(30000), A(-30000)) == 5536);
    CHECK(angle16_distance(A(0), A(-32768)) == 32768);
    CHECK(angle16_distance(A(5), A(5)) == 0);

    for (i = -32768; i <= 32767; i += 1021) {
        for (j = -32768; j <= 32767; j += 773) {
            angle16_t a = A(i), b = A(j);
            CHECK(angle16_add(angle16_sub(a, b), b).raw == a.raw);
            CHECK(angle16_add(a, angle16_neg(a)).raw == 0);
            CHECK(angle16_lerp(a, b, 0).raw == a.raw);
            CHECK(angle16_lerp(a, b, 65536).raw == b.raw);
            CHECK(angle16_distance(a, b) == angle16_distance(b, a));
            CHECK(angle16_distance(a, b) == (uint16_t)abs(angle16_sub(a, b).raw)
                  || angle16_sub(a, b).raw == -32768);
        }
    }
    /* lerp end points exactly, for every distance including the wrap. */
    for (i = -32768; i <= 32767; i += 97) {
        angle16_t a = A(-32768), b = A(i);
        CHECK(angle16_lerp(a, b, 65536).raw == b.raw);
        CHECK(angle16_lerp(b, a, 65536).raw == a.raw);
    }
    CHECK(angle16_lerp(A(30000), A(-30000), 32768).raw == -32768);   /* midpoint across the wrap */
    CHECK(angle16_lerp(A(0), A(16384), 32768).raw == 8192);
    CHECK(angle16_lerp(A(10), A(-10), 32768).raw == 0);
    CHECK(angle16_lerp(A(-10), A(10), 32768).raw == 0);
}

/* sin(u * pi/32768) with the quarter-wave fold done exactly, so that the
 * reference is exactly 0 at pi rather than sin(M_PI) = 1.2e-16. */
static double sinpi_ref(int r)
{
    unsigned k = angle16__quarter((uint16_t)r);
    double s = sin(k * PI_D / 32768.0);
    return ((uint16_t)r & 32768u) ? -s : s;
}

/* Error of got against the true value, in ulps of the float nearest true. */
static double ulp_error(float got, double truth)
{
    double mag = fabs(truth), ulp;
    if (mag == 0.0) return got == 0.0f ? 0.0 : 1e9;
    ulp = ldexp(1.0, (int)floor(log2(mag)) - 23);
    return fabs((double)got - truth) / ulp;
}

static void test_trig(void)
{
    double worst = 0.0;
    int r, not_rounded = 0;
    for (r = -32768; r <= 32767; r++) {
        angle16_t a = A(r);
        float s = angle16_sin(a), c = angle16_cos(a);
        double e = ulp_error(s, sinpi_ref(r));
        if (e > worst) worst = e;
        if (e > 0.5) not_rounded++;
        CHECK(e <= 1.25);
        CHECK(same_bits(angle16_sin(angle16_neg(a)), -s));                /* odd  */
        CHECK(same_bits(angle16_cos(angle16_neg(a)), c));                 /* even */
        CHECK(same_bits(angle16_sin(angle16_add(a, ANGLE16_HALF_PI)), c));
        CHECK(same_bits(angle16_sin(angle16_add(a, ANGLE16_PI)), -s));
        CHECK(fabsf(s * s + c * c - 1.0f) < 2.4e-7f);
        CHECK(angle16_atan2(s, c).raw == r);
    }
    /* Every 8-bit angle is a multiple of 256 here, where the fine term
     * vanishes and the result is a coarse table entry: correctly rounded. */
    for (r = -32768; r <= 32767; r += 256)
        CHECK(ulp_error(angle16_sin(A(r)), sinpi_ref(r)) <= 0.5);
    CHECK(angle16_sin(A(0)) == 0.0f);
    CHECK(angle16_sin(A(16384)) == 1.0f);
    CHECK(angle16_sin(A(-16384)) == -1.0f);
    CHECK(angle16_sin(A(-32768)) == 0.0f);      /* exactly, unlike sinf(pi) */
    CHECK(!signbit(angle16_sin(A(-32768))));    /* and it is +0, not -0 */
    CHECK(!signbit(angle16_cos(A(16384))));
    CHECK(angle16_cos(A(0)) == 1.0f);
    CHECK(angle16_cos(A(16384)) == 0.0f);
    CHECK(angle16_cos(A(-16384)) == 0.0f);
    CHECK(angle16_cos(A(-32768)) == -1.0f);
    CHECK(angle16_sin(A(8192)) == angle16_cos(A(8192)));
    printf("  sin worst %.3f ulp, %d of 65536 not correctly rounded\n", worst, not_rounded);

    CHECK(angle16_atan2(1.0f, 0.0f).raw == 16384);
    CHECK(angle16_atan2(0.0f, -1.0f).raw == -32768);
    CHECK(angle16_atan2(0.0f, 1.0f).raw == 0);
    CHECK(angle16_atan2(-1.0f, 0.0f).raw == -16384);
    CHECK(angle16_atan2(1.0f, 1.0f).raw == 8192);
}

static void test_degrees_int(void)
{
    static const int in[]  = { 0, 90, 180, -180, 270, -90, 360, 45, -45, 450, -450,
                               1, 2, 3, 179, -179, 181, 720090 };
    static const int out[] = { 0, 16384, -32768, -32768, -16384, -16384, 0, 8192, -8192, 16384, -16384,
                               182, 364, 546, 32586, -32586, -32586, 16384 };
    size_t i;
    for (i = 0; i < sizeof in / sizeof in[0]; i++)
        CHECK(angle16_from_degrees_int(in[i]).raw == out[i]);
    for (i = 0; i < 3600; i++) {
        int d = (int)i - 1800;
        int ref = angle16_from_degrees((float)d).raw;
        int got = angle16_from_degrees_int(d).raw;
        CHECK(got == ref || (uint16_t)(got - ref) == 1 || (uint16_t)(ref - got) == 1);
    }
}

int main(void)
{
    printf("angle16 tests: TO_FLOAT_BITS=%d FROM_FLOAT_BITS=%d\n",
           ANGLE16_TO_FLOAT_BITS, ANGLE16_FROM_FLOAT_BITS);
    test_to_float();
    test_from_float();
    test_arithmetic();
    test_trig();
    test_degrees_int();
    printf("%d checks, %d failures\n", checks, fails);
    return fails ? 1 : 0;
}
