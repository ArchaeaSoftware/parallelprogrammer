/* Compile-and-run check for the C++ operator overloads. */
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
