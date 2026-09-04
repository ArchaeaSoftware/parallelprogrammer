/*
 * bench.c
 *
 * Micro-benchmark: cast-based vs bit-trick conversions. Build with `make bench`.
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
#define _POSIX_C_SOURCE 199309L
#include "angle8.h"
#include <stdio.h>
#include <time.h>

#define N 4096
#define REPS 20000

static angle8_t in_a[N];
static float    in_f[N];

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

#define BENCH(name, expr) do {                                              \
    double t0 = now();                                                      \
    for (rep = 0; rep < REPS; rep++)                                        \
        for (i = 0; i < N; i++) acc += (expr);                              \
    printf("%-22s %6.2f ns/op  (acc %g)\n", name,                           \
           (now() - t0) / ((double)N * REPS) * 1e9, (double)acc);           \
    acc = 0;                                                                \
} while (0)

int main(void)
{
    int i, rep;
    uint32_t seed = 1;
    float acc = 0;
    for (i = 0; i < N; i++) {
        seed = seed * 1664525u + 1013904223u;
        in_a[i] = angle8_from_uraw((uint8_t)(seed >> 24));
        in_f[i] = (float)(int32_t)seed * (1.0f / 2147483648.0f) * 100.0f;
    }
    BENCH("to_unit_cast",   angle8_to_unit_cast(in_a[i]));
    BENCH("to_unit_bits",   angle8_to_unit_bits(in_a[i]));
    BENCH("from_unit_cast", (float)angle8_from_unit_cast(in_f[i]).raw);
    BENCH("from_unit_bits", (float)angle8_from_unit_bits(in_f[i]).raw);
    BENCH("sin (table)",    angle8_sin(in_a[i]));
    BENCH("sinf (libm)",    sinf(angle8_to_radians(in_a[i])));
    return 0;
}
