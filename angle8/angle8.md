# angle8: an angle in one byte

`angle8.h` is a header-only angle type for C99 and C++11. An angle is a
signed 8-bit integer scaled by pi, so the byte's whole range covers exactly
one turn:

| raw  | radians      | degrees   |
|------|--------------|-----------|
| -128 | -pi          | -180      |
|  -64 | -pi/2        | -90       |
|    0 | 0            | 0         |
|   32 | pi/4         | 45        |
|   64 | pi/2         | 90        |
|  127 | 127/128 * pi | 178.59375 |

Resolution is 1/256 of a turn, about 1.4 degrees. That is coarse for a
compass, but plenty for a heading in a game entity, a sprite rotation, a
direction packed into a network message, or a phase in a signal table.

## Wrap around is the feature

Because 256 raw units make one full turn, ordinary two's-complement overflow
is angular wrap around. Adding 1 to 127, just under +pi, gives -128, which is
-pi: the same direction. So the operations that usually need range checks
need none:

```c
angle8_t heading = angle8_from_degrees(350.0f);   /* raw -7  */
angle8_t target  = angle8_from_radians(0.3f);      /* raw 12  */
angle8_t err     = angle8_sub(target, heading);    /* raw 19: the shortest turn */
```

`angle8_sub(a, b)` is the shortest signed rotation from b to a, always in
[-pi, pi), because that is what a byte subtraction is. `angle8_add`,
`angle8_neg` and `angle8_scale` (an integer multiple of an angle) are byte
arithmetic too. All of it is computed in unsigned, so the compiler need not concern itself with signed overflow.

The same byte subtraction, read as unsigned, can be used to measure distance.
`angle8_ccw_distance(from, to)` is the counter-clockwise difference, 0..255. `angle8_distance(a, b)`
folds that to the absolute separation, 0..128, where 128 is pi.

`angle8_lerp(a, b, t)` interpolates along the shortest arc, with t in
1/256ths of the way from a to b. It takes the signed difference, scales it,
and adds it back, so interpolating from 170 degrees to -170 degrees goes
through 180, not through 0.

## Conversions

"Unit" means multiples of pi: raw / 128, in [-1, 1). Radians, degrees and
turns are each one multiply away from that.

Going from `int8_t` to `float` is cheap because the scale is a power of two. The
constant pi/128 has the same mantissa as pi, so `angle8_to_radians` costs no
precision beyond pi itself, and `angle8_to_degrees` is exact, since 180/128 is
1.40625.

Going from `float` to `int8_t`, `angle8_from_unit` and its radians, degrees and
turns wrappers scale the input, round to nearest even, and wrap. Any finite
input is accepted; 370 degrees becomes 10 degrees. NaN and infinities map to
0. `angle8_from_degrees_int` does the same for integer degrees with no
floating point at all, for targets that have none.

### Two bit tricks

The header carries two alternative conversion paths, each selected by a macro
and each always available under its own name.

`angle8_to_unit_bits` avoids the `int`-to-`float` instruction. It ORs the byte,
with its sign flipped, into the mantissa of 2.0f. That `float` reads as
`2 + (raw + 128)/128`, which is `3 + raw/128`. Subtracting `3.0f` is exact, so the
result is bit-identical to the plain cast under every rounding mode. This implementation is the default, because it delivers performance parity on x86-64 and is expected to be fastest on targets with no (or slow) `int`-to-`float`.

`angle8_from_unit_bits` avoids the `float`-to-`int` conversion. Adding 1.5 * 2^23
pushes the value into the range where a `float`'s ulp is exactly 1, so the FPU
rounds it to an integer in the low mantissa bits. The low byte of that `float`
is the rounded value mod 256, and the mod is the wrap. It matches `lrintf`
under the default round-to-nearest mode and is opt-in because it depends on
that mode. On x86-64, it makes a difference: `lrintf` is an out-of-line libm call unless you pass `-fno-math-errno`, and the trick is about 4x faster.

## Trigonometry

`angle8_sin` looks up a 65-entry quarter-wave table of sin(k * pi/128) for k
in 0..64, each entry correctly rounded to `float`. Bit 6 of the raw byte
mirrors the index within the quarter wave and bit 7 flips the sign, so the
function is exactly odd, sin(pi) is exactly +0 rather than the 1.2e-16 you
get from `sinf(M_PI)`, and sin(pi/2) is exactly 1. `angle8_cos` is sin of the
angle plus a quarter turn, so it is exactly even and bit-identical to the
sin it is defined by.

Why not convert to radians and call `sinf`? Speed is the obvious reason, but
not the biggest. Measured over all 256 angles on x86-64 with glibc 2.31:

| | `angle8_sin` | `sinf(angle8_to_radians(a))` |
|---|---|---|
| ns per call | 0.55 | 2.5 |
| correctly rounded | 256 of 256 | 147 of 256 |
| worst error where sin is nonzero | 0.5 ulp | 50.5 ulp |
| sin(pi), cos(pi/2) | exactly 0 | 8.74e-8 |
| sin(a + pi) == -sin(a), bit-exact | always | fails for 172 of 256 |
| cos(-a) == cos(a), bit-exact | always | fails for 174 of 256 |

The 50 ulp is not `sinf`'s fault. pi as a `float` is off by 8.7e-8, and near a
half turn that error is comparable to the sine itself, so an angle of 127/128
pi comes back with a relative error of 3e-6. The table never forms the radian
value, so it never inherits pi's rounding. Computing `sinf(x * pi)` from the
unit value has exactly the same problem. A real `sinpif`, which C23 specifies
and glibc has had since 2.41, reduces the argument exactly and would recover
the exact zeros and the symmetries, leaving the table ahead on speed, correct
rounding, and having no libm dependency.

`angle8_atan2(y, x)` is libm's `atan2f` rounded to the nearest raw value. It
inverts sin and cos for all 256 angles.

## C++

The type is a plain struct with one `int8_t` member, `raw`, so it is trivially
copyable and one byte in an array. In C++, the header adds `+ - * += -= *= ==
!=` on top of the C functions, so the usage example above becomes
`target - heading`. `ANGLE8_LIT(r)` makes a constant from a raw value in
either language, and `ANGLE8_ZERO`, `ANGLE8_QUARTER_PI`, `ANGLE8_HALF_PI`,
`ANGLE8_PI` and `ANGLE8_NEG_HALF_PI` are provided.

## Build options

Define before including:

| macro | default | effect |
|---|---|---|
| `ANGLE8_TO_FLOAT_BITS` | 1 | mantissa trick for `int8_t`-to-`float` |
| `ANGLE8_FROM_FLOAT_BITS` | 0 | magic-number trick for `float`-to-`int8_t` |
| `ANGLE8_INLINE` | `static inline` | linkage of every function |

The header needs `<math.h>`, so some C toolchains want `-lm`.

## Testing

The input space is 256 values, so most of the tests are exhaustive rather
than sampled. Every angle round-trips through each conversion, the `float`
conversions are swept across their whole input range including the wrap and
the NaN and infinity cases, a strided grid of angle pairs is added,
subtracted and measured, and every angle's sin and cos are checked against
double precision and against each other for the exact symmetries above. The suite runs under `-Wall -Wextra -pedantic -Werror` in
C99 and C++11 with both trick configurations, 2.4 million checks in all.

## Performance (x86-64, gcc 9, -O2)

| path | ns/op |
|---|---|
| `angle8_to_unit`, either variant | 0.55 |
| `angle8_from_unit_cast` | 2.19 |
| `angle8_from_unit_cast` with `-fno-math-errno` | 0.55 |
| `angle8_from_unit_bits` | 0.61 |
| `angle8_sin` | 0.55 |
| `sinf` on the same angle | 2.5 |

## License

BSD 3-Clause. Copyright (C) 2026 by Archaea Software, LLC.
