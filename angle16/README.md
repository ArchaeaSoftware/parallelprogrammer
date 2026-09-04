# angle16

Header-only 16-bit angle type for C99 and C++11: the 16-bit edition of
[angle8](../angle8/). One `int16_t` holds an angle as a signed integer scaled
by pi:

| raw    | radians           | degrees          |
|--------|-------------------|------------------|
| -32768 | -pi               | -180             |
| -16384 | -pi/2             | -90              |
|      0 | 0                 | 0                |
|  16384 | pi/2              | 90               |
|  32767 | 32767/32768 * pi  | 179.9945068359375 |

Resolution is 1/65536 of a turn, 0.0055 degrees. Two's-complement wrap around
is angular wrap around, so add, subtract, negate and shortest-signed-difference
are plain 16-bit arithmetic with no range checks.

The API is angle8's with the prefix changed: `angle16_from_degrees`,
`angle16_sub`, `angle16_sin` and so on, with the same C++ operators. This file
covers what differs; see the angle8 README for the rest.

## What changes at 16 bits

Most of the header is a width substitution. The properties that make the 8-bit
design pleasant survive:

- Degrees stay exact: 180/32768 is 45 * 2^-13, and raw * 45 fits in a `float`
  mantissa.
- The `int16_t`-to-`float` mantissa trick works unchanged with the sixteen bits
  placed at mantissa bit 7 instead of bit 15, and is still bit-identical to
  the cast.
- The `float`-to-`int16_t` magic-number trick works unchanged, but its range is
  128 turns rather than 32768, because the value has to stay below 2^22
  before the 1.5 * 2^23 is added. Past that it falls back to the conversion.

Three things had to be different:

- `angle16_from_unit_cast` uses `llrintf`, and the "so large it must be a
  whole number of turns" cutoff moves from 2^31 to 2^39, where a `float`'s ulp
  reaches 65536.
- `angle16_lerp` takes t in 1/65536ths, and the product of the signed
  difference and t reaches 2^31 in magnitude, so it is formed in 64 bits.
- sin cannot be a quarter-wave table: 16385 floats is 64 KB.

## sin

After the quarter-wave fold, sin only ever sees k * pi/32768 with k in
0..16384. The index is split as k = 128 hi + lo and the angle-addition
identity is evaluated as

    sin(hi + lo) = sin(hi) + (cos(hi) sin(lo) - sin(hi) (1 - cos(lo)))

from a 129-entry coarse table of sin(hi * pi/256), which also serves cos(hi)
read from the other end, and two 128-entry fine tables of sin(lo * pi/32768)
and 1 - cos(lo * pi/32768). That is 385 floats, 1540 bytes. Storing 1 - cos
rather than cos keeps the whole correction small, at most 1.3 percent of the
result, so the answer is a correctly rounded coarse entry plus a small term
rounded once.

Measured over all 65536 angles:

| | |
|---|---|
| worst error | 1.21 ulp |
| correctly rounded | 49332 of 65536 (75 percent) |
| over 1 ulp | 24 of 65536, all at coarse index 1, 4 or 13 |
| multiples of pi/128, the angles angle8 has | correctly rounded, the fine term vanishes |
| sin(pi), cos(pi/2) | exactly +0 |
| sin(pi/2) | exactly 1 |
| sin odd, cos even, sin(a + pi/2) == cos(a) | bit-exact |

The worst cases sit where the coarse entry and the correction are the same
size, so three roundings stack. A double-precision coarse table would fix
that at the cost of double arithmetic, which the targets that want the bit
tricks tend not to have. As it stands the accuracy is the same class as a
good libm `sinf`, and unlike `sinf` on a converted radian value it never
inherits the rounding of pi.

## Build and test

```
make test    # C99 and C++11, -Wall -Wextra -pedantic -Werror, both trick configs
make bench   # ns/op for each conversion path and sin
```

The tests are exhaustive over all 65536 angles for the conversions, the
round trips and sin, with the arithmetic checked on a strided grid of pairs.

## Performance (x86-64, gcc 9, -O2)

| path | ns/op |
|---|---|
| `to_unit_cast` | 0.63 |
| `to_unit_bits` | 0.55 |
| `from_unit_cast` | 2.19 |
| `from_unit_bits` | 0.61 |
| `angle16_sin` | 0.78 |
| `sinf` on the same angle | 2.4 |

## License

BSD 3-Clause. Copyright (C) 2026 by Archaea Software, LLC. See [LICENSE](LICENSE).
