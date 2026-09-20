# Expected benchmark output

This directory holds a reference run of the drivers in `bench/`, produced by
`bench/replicate.sh <build-dir> bench/expected` on an idle machine.

It is deliberately empty in this checkout. An earlier run was discarded
because the machine was not quiet: another CUDA job held the GPU at 100%
utilization and two processes were consuming CPU cores, which roughly halved
every figure. `replicate.sh` now refuses to run under those conditions.

To populate it, build the project and run:

    bench/replicate.sh build bench/expected

on a machine with no other GPU or CPU load, with the GPU at its default power
limit. `machine.txt` records the hardware, the load average, the compiler, and
the GPU's power limit and maximum clocks, all of which the timings depend on.

Timings are hardware-dependent and will not match elsewhere. What a
replication should confirm is the structure, not the digits:

- appending a limb-column costs about the same whatever the column's current
  width, while a rescale costs more as the column grows, because it rewrites
  every limb-column the column holds;
- the CPU scales with threads while the working set fits in cache, and stops
  scaling once the input exceeds last-level cache;
- batching several matrices into one pass helps most where accumulation matrix
  traffic is the bound;
- device-resident input beats input streamed from the host, which is limited
  by the bus.
