/*
 *
 * asan_test.cpp
 *
 * Correctness and bounds check for the masked tails.
 *
 * Build with: g++ -std=c++17 -O1 -g -fsanitize=address -mavx512f -mavx512dq -mavx512bw -mavx512vl -mfma asan_test.cpp
 *
 * Copyright (C) 2026 by Nicholas Wilt.
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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND
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

// Correctness and bounds check for the masked tails.
//
// Buffers are exactly n floats, so AddressSanitizer will flag any access past
// the end. Sizes deliberately include n < 16, n == 16, and awkward remainders.
#include "ops.hpp"
#include "fusion.hpp"

#include <cstdio>
#include <vector>

int main() {
    const std::size_t sizes[] = {1, 2, 15, 16, 17, 31, 33, 63, 65,
                                 127, 129, 1000, 4095, 65537};
    float coeffs[16];
    for (int i = 0; i < 16; ++i) coeffs[i] = 0.1f * float(i + 1);

    int bad = 0;
    for (std::size_t n : sizes) {
        std::vector<float> in(n), out(n), tmp(n);
        for (std::size_t i = 0; i < n; ++i)
            in[i] = 0.5f + 0.001f * float(i);

        Scale   s(0.5f);
        FmaBias f(1.25f, 0.125f);
        Clamp   cl(-8.0f, 8.0f);
        Poly<5> p(coeffs);
        Relu    r(0.0f);

        auto sum = [&] {
            double a = 0;
            for (std::size_t i = 0; i < n; ++i) a += out[i];
            return a;
        };

        separate_passes(in.data(), out.data(), tmp.data(), n, s, f, cl, p, r);
        double a = sum();
        fused_vecmajor<4>(in.data(), out.data(), n, s, f, cl, p, r);
        double b = sum();
        fused_vecmajor_seq<4>(in.data(), out.data(), n, s, f, cl, p, r);
        double c = sum();
        fused_opmajor_byval<4>(in.data(), out.data(), n, s, f, cl, p, r);
        double d = sum();
        fused_opmajor<4>(in.data(), out.data(), n, s, f, cl, p, r);
        double e = sum();

        bool ok = (a == b) && (b == c) && (c == d) && (d == e);
        if (!ok) ++bad;
        std::printf("n=%6zu  %s\n", n, ok ? "ok" : "MISMATCH");
    }
    std::printf("%s\n", bad ? "FAILURES" : "all sizes agree");
    return bad != 0;
}
