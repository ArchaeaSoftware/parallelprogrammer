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
#include "angle8.h"
#include <cstdio>

int main()
{
    int fails = 0;
    angle8_t a = ANGLE8_LIT(120);
    angle8_t b = ANGLE8_LIT(-120);

    if ((a - b) != ANGLE8_LIT(-16)) fails++;
    if ((b - a) != ANGLE8_LIT(16))  fails++;
    if (-ANGLE8_PI != ANGLE8_PI)    fails++;
    if (ANGLE8_HALF_PI * 2 != ANGLE8_PI) fails++;
    if (2 * ANGLE8_HALF_PI != ANGLE8_PI) fails++;

    angle8_t c = ANGLE8_LIT(127);
    c += ANGLE8_LIT(1);
    if (c != ANGLE8_PI) fails++;
    c -= ANGLE8_LIT(1);
    if (c != ANGLE8_LIT(127)) fails++;
    c *= 2;
    if (c != ANGLE8_LIT(-2)) fails++;

    if (angle8_to_degrees(ANGLE8_HALF_PI) != 90.0f) fails++;
    if (angle8_from_degrees(450.0f) != ANGLE8_HALF_PI) fails++;

    std::printf("angle8 C++ tests: %d failures\n", fails);
    return fails ? 1 : 0;
}
