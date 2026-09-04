/*
 * test_cpp.cpp
 *
 * Compile-and-run check for the C++ operator overloads.
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
#include <cstdio>

int main()
{
    int fails = 0;
    angle16_t a = ANGLE16_LIT(30000);
    angle16_t b = ANGLE16_LIT(-30000);

    if ((a - b) != ANGLE16_LIT(-5536)) fails++;
    if ((b - a) != ANGLE16_LIT(5536))  fails++;
    if (-ANGLE16_PI != ANGLE16_PI)     fails++;
    if (ANGLE16_HALF_PI * 2 != ANGLE16_PI) fails++;
    if (2 * ANGLE16_HALF_PI != ANGLE16_PI) fails++;

    angle16_t c = ANGLE16_LIT(32767);
    c += ANGLE16_LIT(1);
    if (c != ANGLE16_PI) fails++;
    c -= ANGLE16_LIT(1);
    if (c != ANGLE16_LIT(32767)) fails++;
    c *= 2;
    if (c != ANGLE16_LIT(-2)) fails++;

    if (angle16_to_degrees(ANGLE16_HALF_PI) != 90.0f) fails++;
    if (angle16_from_degrees(450.0f) != ANGLE16_HALF_PI) fails++;

    std::printf("angle16 C++ tests: %d failures\n", fails);
    return fails ? 1 : 0;
}
