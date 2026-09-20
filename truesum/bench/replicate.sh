#!/bin/sh
# Runs every benchmark driver and writes its output next to the expected
# output shipped with the software, so a reviewer can diff the two. Timings
# are hardware-dependent; the expected files record the machine they came
# from in their first line.
#
# Usage: bench/replicate.sh <build-dir> [output-dir]
set -e
build=${1:?usage: replicate.sh <build-dir> [output-dir]}
out=${2:-bench/output}
mkdir -p "$out"

# Record the machine, since every timing below is hardware-dependent. The GPU
# power limit and clocks are part of that record: a capped board streams
# proportionally less, and every kernel here is bandwidth-bound.
{
    echo "# replicated $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "# cpu:    $(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2- | sed 's/^ *//')"
    echo "# cores:  $(nproc 2>/dev/null) logical"
    echo "# load:   $(cut -d' ' -f1-3 /proc/loadavg 2>/dev/null)"
    echo "# os:     $(uname -sr)"
    echo "# cxx:    $(${CXX:-c++} --version 2>/dev/null | head -1)"
    if command -v nvidia-smi >/dev/null 2>&1; then
        echo "# gpu:    $(nvidia-smi --query-gpu=name,driver_version --format=csv,noheader 2>/dev/null | head -1)"
        echo "# power:  $(nvidia-smi --query-gpu=power.limit,power.default_limit --format=csv,noheader 2>/dev/null | head -1) (limit, default)"
        echo "# clocks: $(nvidia-smi --query-gpu=clocks.max.sm,clocks.max.mem --format=csv,noheader 2>/dev/null | head -1) (max sm, max mem)"
    fi
    command -v nvcc >/dev/null 2>&1 && \
        echo "# nvcc:   $(nvcc --version | grep release | sed 's/^ *//')"
} > "$out/machine.txt"
cat "$out/machine.txt"

# Refuse to publish numbers measured on a contended machine. Every driver here
# is bandwidth-bound, so a co-resident job does not add a few percent of noise
# -- it halves the result. Set TRUESUM_BENCH_FORCE=1 to measure anyway.
contended=""
load=$(cut -d' ' -f1 /proc/loadavg 2>/dev/null | cut -d. -f1)
[ -n "$load" ] && [ "$load" -ge 2 ] && \
    contended="$contended\n  load average is $(cut -d' ' -f1 /proc/loadavg); another job is using the CPU"
if command -v nvidia-smi >/dev/null 2>&1; then
    others=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null | grep -c . || true)
    [ "${others:-0}" -gt 0 ] && \
        contended="$contended\n  $others other process(es) hold GPU memory; the card is not idle"
    util=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits 2>/dev/null | head -1)
    [ -n "$util" ] && [ "$util" -ge 20 ] && \
        contended="$contended\n  GPU utilization is ${util}% before any benchmark started"
    lim=$(nvidia-smi --query-gpu=power.limit --format=csv,noheader,nounits 2>/dev/null | head -1 | cut -d. -f1)
    def=$(nvidia-smi --query-gpu=power.default_limit --format=csv,noheader,nounits 2>/dev/null | head -1 | cut -d. -f1)
    [ -n "$lim" ] && [ -n "$def" ] && [ "$lim" -lt "$def" ] && \
        contended="$contended\n  GPU power limit is ${lim} W against a ${def} W default; clocks will be lower"
fi
if [ -n "$contended" ]; then
    printf '\nWILL NOT MEASURE. The machine is not quiet:%b\n\n' "$contended"
    if [ "${TRUESUM_BENCH_FORCE:-0}" != "1" ]; then
        echo "Re-run when the machine is idle, or set TRUESUM_BENCH_FORCE=1 to override."
        exit 1
    fi
    echo "TRUESUM_BENCH_FORCE=1: measuring anyway. These numbers are not comparable."
fi

run() {
    name=$1; shift
    bin=${name%_scalar_kernel}
    if [ -x "$build/$bin" ]; then
        echo "== $name"
        "$build/$bin" "$@" | tee "$out/$name.txt"
    else
        echo "== $name: not built (skipped)"
    fi
}

run bench_table_cpu
run bench_adjust_cpu
run bench_kernels
TRUESUM_KERNEL=scalar run bench_kernels_scalar_kernel
run bench_table_gpu
run bench_adjust_gpu

echo
echo "outputs in $out/; expected outputs in bench/expected/"
