# angle8

Header-only 8-bit angle type for C99 and C++11. One byte holds an angle as a
signed integer scaled by pi:

| raw  | radians        | degrees   |
|------|----------------|-----------|
| -128 | -pi            | -180      |
|  -64 | -pi/2          | -90       |
|    0 | 0              | 0         |
|   64 | pi/2           | 90        |
|  127 | 127/128 * pi   | 178.59375 |

Two's-complement wrap around is angular wrap around, so add, subtract, negate
and shortest-signed-difference are plain byte arithmetic with no range checks.

## Usage

```c
#include "angle8.h"

angle8_t heading = angle8_from_degrees(350.0f);       /* raw -7 */
angle8_t target  = angle8_from_radians(0.3f);          /* raw 12 */
angle8_t err     = angle8_sub(target, heading);        /* raw 19: shortest turn */
float    rad     = angle8_to_radians(err);
float    s       = angle8_sin(heading);                /* table lookup */
```

In C++ the struct also has `+ - * += -= *= == !=` overloads.

## API

| Function | Notes |
|---|---|
| `angle8_from_raw / angle8_raw / angle8_uraw / angle8_from_uraw` | raw int8 / uint8 access |
| `angle8_to_unit / _radians / _degrees / _turns` | one multiply each; degrees and turns are exact |
| `angle8_from_unit / _radians / _degrees / _turns` | wraps any finite input, rounds to nearest even; NaN and inf give 0 |
| `angle8_from_degrees_int` | integer only, no FPU |
| `angle8_add / _sub / _neg / _scale / _eq` | mod 2pi |
| `angle8_ccw_distance / angle8_distance` | 0..255 counter-clockwise, or 0..128 absolute |
| `angle8_lerp(a, b, t)` | shortest arc, `t` in 1/256ths |
| `angle8_sin / _cos` | 65-entry quarter-wave float table, exact symmetry |
| `angle8_atan2` | libm atan2f rounded to the nearest raw value |

"Unit" means multiples of pi: raw / 128, in [-1, 1).

## Float conversion tricks

`angle8_to_unit_bits` ORs the byte into the mantissa of 2.0f and subtracts
3.0f. The result is bit-identical to the plain cast under every rounding
mode. It is the default (`ANGLE8_TO_FLOAT_BITS=1`).

`angle8_from_unit_bits` adds 1.5 * 2^23 so the FPU rounds to an integer in the
low mantissa bits, then reads the low byte, which also performs the mod-256
wrap. It matches `lrintf` under the default round-to-nearest mode. It is
opt-in (`ANGLE8_FROM_FLOAT_BITS=1`) because it depends on the rounding mode.

Both variants are always available under their explicit names so you can pick
per call site or benchmark them with `make bench`.

## Build and test

```
make test    # C99 and C++11, -Wall -Wextra -pedantic -Werror, both trick configs
make bench   # ns/op for each conversion path
```

## Performance notes (x86-64, gcc 9 / clang 10, -O2)

| path | ns/op | comment |
|---|---|---|
| `to_unit_cast` | 0.55 | `cvtsi2ss` + `mulss` |
| `to_unit_bits` | 0.55 | `add, movzx, shl, or, movd, subss`; no int-to-float unit needed |
| `from_unit_cast` | 2.19 | `lrintf` is an out-of-line libm call unless you pass `-fno-math-errno` |
| `from_unit_cast` with `-fno-math-errno` | 0.55 | inlines to `cvtss2si` |
| `from_unit_bits` | 0.61 | `addss` + read low byte, no libm |
| `angle8_sin` (table) | 0.55 | vs 2.5 for `sinf` |

So on x86-64 the int8-to-float trick is a wash; it pays off on targets without a
fast int-to-float instruction. The float-to-int8 trick is the easy win when
you cannot or do not want to pass `-fno-math-errno`.

## License

BSD 3-Clause. Copyright (C) 2026 by Archaea Software, LLC. See [LICENSE](LICENSE).
