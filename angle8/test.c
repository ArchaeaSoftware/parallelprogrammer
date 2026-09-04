/*
 * test.c
 *
 * Exhaustive tests for angle8.h. Build and run with `make test`.
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
#include "angle8.h"
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
static angle8_t A(int r) { return angle8_from_raw((int8_t)r); }

static void test_to_float(void)
{
    int r;
    for (r = -128; r <= 127; r++) {
        angle8_t a = A(r);
        float unit = angle8_to_unit(a);
        double ref = r * PI_D / 128.0;
        float rad = angle8_to_radians(a);
        CHECK(same_bits(angle8_to_unit_bits(a), angle8_to_unit_cast(a)));
        CHECK(unit == (float)r / 128.0f);
        CHECK(fabs(rad - ref) <= fabs(ref) * 1.2e-7);   /* within 1 float ulp */
        CHECK(angle8_to_degrees(a) == (float)r * 1.40625f);
        CHECK(angle8_to_turns(a) == (float)r / 256.0f);
    }
    /* Exact endpoints and landmarks. */
    CHECK(angle8_to_radians(A(-128)) == -ANGLE8_PI_F);
    CHECK(angle8_to_radians(A(127)) == (127.0f / 128.0f) * ANGLE8_PI_F);
    CHECK(angle8_to_radians(A(64)) == ANGLE8_PI_F * 0.5f);
    CHECK(angle8_to_radians(A(0)) == 0.0f);
    CHECK(angle8_to_degrees(A(64)) == 90.0f);
    CHECK(angle8_to_degrees(A(-64)) == -90.0f);
    CHECK(angle8_to_degrees(A(-128)) == -180.0f);
    CHECK(angle8_to_degrees(A(127)) == 178.59375f);
    CHECK(angle8_to_degrees(A(1)) == 1.40625f);
    CHECK(angle8_to_turns(A(64)) == 0.25f);
    CHECK(angle8_to_turns(A(-128)) == -0.5f);
}

static void test_from_float(void)
{
    int r, j, k;
    uint32_t seed = 12345;

    /* Round trips through every representation, both implementations. */
    for (r = -128; r <= 127; r++) {
        angle8_t a = A(r);
        CHECK(angle8_from_radians(angle8_to_radians(a)).raw == r);
        CHECK(angle8_from_degrees(angle8_to_degrees(a)).raw == r);
        CHECK(angle8_from_turns(angle8_to_turns(a)).raw == r);
        CHECK(angle8_from_unit_cast(angle8_to_unit(a)).raw == r);
        CHECK(angle8_from_unit_bits(angle8_to_unit(a)).raw == r);
    }

    /* Bit-trick and cast paths agree exactly on exact ties (k/256 -> x = k/2). */
    for (k = -200000; k <= 200000; k++) {
        float u = (float)k / 256.0f;
        CHECK(angle8_from_unit_bits(u).raw == angle8_from_unit_cast(u).raw);
    }
    /* ...and on pseudo-random values across the trick's full range. */
    for (j = 0; j < 2000000; j++) {
        float u;
        seed = seed * 1664525u + 1013904223u;
        memcpy(&u, &seed, 4);
        if (!(fabsf(u) < 40000.0f)) u = (float)(seed >> 8) * (1.0f / 16777216.0f) * 8192.0f - 4096.0f;
        CHECK(angle8_from_unit_bits(u).raw == angle8_from_unit_cast(u).raw);
    }

    /* Ties round to even. */
    CHECK(angle8_from_unit( 0.5f / 128.0f).raw == 0);
    CHECK(angle8_from_unit( 1.5f / 128.0f).raw == 2);
    CHECK(angle8_from_unit(-0.5f / 128.0f).raw == 0);
    CHECK(angle8_from_unit(-1.5f / 128.0f).raw == -2);

    /* Wrap: adding whole turns changes nothing (sums chosen to be exact). */
    for (j = -512; j <= 512; j++) {
        float t = (float)j / 1024.0f;
        angle8_t base = angle8_from_turns(t);
        for (k = -3; k <= 3; k++) {
            CHECK(angle8_from_turns(t + (float)k).raw == base.raw);
            CHECK(angle8_from_unit_bits((t + (float)k) * 2.0f).raw == base.raw);
        }
    }
    CHECK(angle8_from_unit(1.0f).raw == -128);
    CHECK(angle8_from_unit(-1.0f).raw == -128);
    CHECK(angle8_from_unit(2.0f).raw == 0);
    CHECK(angle8_from_unit(3.0f).raw == -128);
    CHECK(angle8_from_degrees(360.0f).raw == 0);
    CHECK(angle8_from_degrees(270.0f).raw == -64);
    CHECK(angle8_from_degrees(-90.0f).raw == -64);
    CHECK(angle8_from_radians(ANGLE8_PI_F).raw == -128);
    CHECK(angle8_from_radians(-ANGLE8_PI_F).raw == -128);

    /* Non-finite and huge inputs are well defined. */
    CHECK(angle8_from_radians(NAN).raw == 0);
    CHECK(angle8_from_radians(INFINITY).raw == 0);
    CHECK(angle8_from_radians(-INFINITY).raw == 0);
    CHECK(angle8_from_unit_bits(NAN).raw == 0);
    CHECK(angle8_from_unit(1e30f).raw == 0);
    CHECK(angle8_from_unit(1e8f).raw == 0);                       /* |x| >= 2^31: multiple of 256 */
    CHECK(angle8_from_unit(40000.0f).raw == 0);                   /* x = 5120000 = 20000 * 256 */
    CHECK(angle8_from_unit(40000.0078125f).raw == 1);             /* x = 5120001 */
    CHECK(angle8_from_unit_bits(40000.0078125f).raw == 1);        /* falls back past 2^22 */
    CHECK(angle8_from_unit_bits(1e7f).raw == angle8_from_unit_cast(1e7f).raw);
}

static void test_arithmetic(void)
{
    int i, j;
    CHECK(angle8_add(A(127), A(1)).raw == -128);
    CHECK(angle8_add(A(-128), A(-1)).raw == 127);
    CHECK(angle8_sub(A(-128), A(1)).raw == 127);
    CHECK(angle8_sub(A(127), A(-1)).raw == -128);
    CHECK(angle8_neg(A(-128)).raw == -128);
    CHECK(angle8_neg(A(0)).raw == 0);
    CHECK(angle8_neg(A(1)).raw == -1);
    CHECK(angle8_scale(A(64), 2).raw == -128);
    CHECK(angle8_scale(A(64), 4).raw == 0);
    CHECK(angle8_scale(A(1), 256).raw == 0);
    CHECK(angle8_scale(A(3), -1).raw == -3);
    CHECK(angle8_scale(A(-128), -1).raw == -128);
    CHECK(angle8_eq(ANGLE8_PI, A(-128)));
    CHECK(ANGLE8_HALF_PI.raw == 64 && ANGLE8_NEG_HALF_PI.raw == -64);
    CHECK(ANGLE8_QUARTER_PI.raw == 32 && ANGLE8_ZERO.raw == 0);

    /* Shortest signed difference crosses the wrap. */
    CHECK(angle8_sub(A(-120), A(120)).raw == 16);
    CHECK(angle8_sub(A(120), A(-120)).raw == -16);
    CHECK(angle8_ccw_distance(A(120), A(-120)) == 16);
    CHECK(angle8_ccw_distance(A(-120), A(120)) == 240);
    CHECK(angle8_distance(A(120), A(-120)) == 16);
    CHECK(angle8_distance(A(0), A(-128)) == 128);
    CHECK(angle8_distance(A(5), A(5)) == 0);

    for (i = -128; i <= 127; i += 7) {
        for (j = -128; j <= 127; j += 5) {
            angle8_t a = A(i), b = A(j);
            CHECK(angle8_add(angle8_sub(a, b), b).raw == a.raw);
            CHECK(angle8_add(a, angle8_neg(a)).raw == 0);
            CHECK(angle8_lerp(a, b, 0).raw == a.raw);
            CHECK(angle8_lerp(a, b, 256).raw == b.raw);
            CHECK(angle8_distance(a, b) == angle8_distance(b, a));
            CHECK(angle8_distance(a, b) == (uint8_t)abs(angle8_sub(a, b).raw)
                  || angle8_sub(a, b).raw == -128);
        }
    }
    CHECK(angle8_lerp(A(120), A(-120), 128).raw == -128);   /* midpoint across the wrap */
    CHECK(angle8_lerp(A(0), A(64), 128).raw == 32);
    CHECK(angle8_lerp(A(10), A(-10), 128).raw == 0);
    CHECK(angle8_lerp(A(-10), A(10), 128).raw == 0);
}

static void test_trig(void)
{
    int r;
    for (r = -128; r <= 127; r++) {
        angle8_t a = A(r);
        double ang = r * PI_D / 128.0;
        float s = angle8_sin(a), c = angle8_cos(a);
        CHECK(fabs(s - sin(ang)) <= 1.2e-7);
        CHECK(fabs(c - cos(ang)) <= 1.2e-7);
        CHECK(same_bits(angle8_sin(angle8_neg(a)), -s));                 /* odd  */
        CHECK(same_bits(angle8_cos(angle8_neg(a)), c));                  /* even */
        CHECK(same_bits(angle8_sin(angle8_add(a, ANGLE8_HALF_PI)), c));
        CHECK(same_bits(angle8_sin(angle8_add(a, ANGLE8_PI)), -s));
        CHECK(fabsf(s * s + c * c - 1.0f) < 2.4e-7f);
        CHECK(angle8_atan2(s, c).raw == r);
    }
    CHECK(angle8_sin(A(0)) == 0.0f);
    CHECK(angle8_sin(A(64)) == 1.0f);
    CHECK(angle8_sin(A(-64)) == -1.0f);
    CHECK(angle8_sin(A(-128)) == 0.0f);      /* exactly, unlike sinf(pi) */
    CHECK(!signbit(angle8_sin(A(-128))));    /* and it is +0, not -0 */
    CHECK(!signbit(angle8_cos(A(64))));
    CHECK(angle8_cos(A(0)) == 1.0f);
    CHECK(angle8_cos(A(64)) == 0.0f);
    CHECK(angle8_cos(A(-64)) == 0.0f);
    CHECK(angle8_cos(A(-128)) == -1.0f);
    CHECK(angle8_sin(A(32)) == angle8_cos(A(32)));

    CHECK(angle8_atan2(1.0f, 0.0f).raw == 64);
    CHECK(angle8_atan2(0.0f, -1.0f).raw == -128);
    CHECK(angle8_atan2(0.0f, 1.0f).raw == 0);
    CHECK(angle8_atan2(-1.0f, 0.0f).raw == -64);
    CHECK(angle8_atan2(1.0f, 1.0f).raw == 32);
}

static void test_degrees_int(void)
{
    static const int in[]  = { 0, 90, 180, -180, 270, -90, 360, 45, -45, 450, -450,
                               1, 2, 3, 179, -179, 181, 720090 };
    static const int out[] = { 0, 64, -128, -128, -64, -64, 0, 32, -32, 64, -64,
                               1, 1, 2, 127, -127, -127, 64 };
    size_t i;
    for (i = 0; i < sizeof in / sizeof in[0]; i++)
        CHECK(angle8_from_degrees_int(in[i]).raw == out[i]);
    for (i = 0; i < 3600; i++) {
        int d = (int)i - 1800;
        int ref = angle8_from_degrees((float)d).raw;
        int got = angle8_from_degrees_int(d).raw;
        CHECK(got == ref || (uint8_t)(got - ref) == 1 || (uint8_t)(ref - got) == 1);
    }
}

int main(void)
{
    printf("angle8 tests: TO_FLOAT_BITS=%d FROM_FLOAT_BITS=%d\n",
           ANGLE8_TO_FLOAT_BITS, ANGLE8_FROM_FLOAT_BITS);
    test_to_float();
    test_from_float();
    test_arithmetic();
    test_trig();
    test_degrees_int();
    printf("%d checks, %d failures\n", checks, fails);
    return fails ? 1 : 0;
}
