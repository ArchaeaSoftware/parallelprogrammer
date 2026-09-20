#!/bin/sh
# Assembles the software component of a TOMS Algorithm submission in the
# layout the TOMS Algorithms guidelines (v1.0, November 2024) ask for:
#
#   userManual.pdf
#   userManual/        LaTeX source of the manual
#   CPP/               the C++ library, tests, demo and CPU benchmarks
#   CUDA/              the CUDA class, its test and GPU benchmarks
#   Doc/               LICENSE, README and top-level build files
#
# No .git, no build products. Run from the truesum directory:
#   tools/make_calgo_archive.sh [output.zip]
set -e
cd "$(dirname "$0")/.."
out=${1:-truesum-calgo.zip}
stage=$(mktemp -d)
trap 'rm -rf "$stage"' EXIT

if [ ! -f userManual/userManual.pdf ]; then
    echo "userManual/userManual.pdf is missing; build it first:" >&2
    echo "  (cd userManual && latexmk -pdf userManual.tex)" >&2
    exit 1
fi
cp userManual/userManual.pdf "$stage/"
mkdir -p "$stage/userManual" "$stage/CPP" "$stage/CUDA" "$stage/Doc"
cp userManual/userManual.tex "$stage/userManual/"

# CPP/: everything the CPU build needs, in the tree layout CMakeLists expects.
for d in include src tests examples bench; do
    mkdir -p "$stage/CPP/$d"
done
cp include/truesum/*.hpp "$stage/CPP/include/truesum/" 2>/dev/null || {
    mkdir -p "$stage/CPP/include/truesum" && cp include/truesum/*.hpp "$stage/CPP/include/truesum/"; }
cp src/*.cpp src/*.hpp "$stage/CPP/src/"
cp tests/test_truesum.cpp "$stage/CPP/tests/"
cp examples/demo.cpp "$stage/CPP/examples/"
cp bench/bench_common.hpp bench/*.cpp bench/replicate.sh "$stage/CPP/bench/"
cp -r bench/expected "$stage/CPP/bench/" 2>/dev/null || true
cp CMakeLists.txt "$stage/CPP/"

# CUDA/: the device class and what exercises it.
mkdir -p "$stage/CUDA/src" "$stage/CUDA/tests" "$stage/CUDA/bench"
cp src/cuda_accumulation_matrix.cu "$stage/CUDA/src/"
cp include/truesum/cuda_accumulation_matrix.hpp "$stage/CUDA/"
cp tests/test_cuda.cu "$stage/CUDA/tests/"
cp bench/*.cu "$stage/CUDA/bench/"

# Doc/: license and readme. The BSD-3 license text lives at the repository
# root, one level above this directory.
cp ../LICENSE "$stage/Doc/LICENSE"
cp README.md "$stage/Doc/README.md"
cp docs/*.md "$stage/Doc/"

rm -f "$out"
(cd "$stage" && zip -qr "$OLDPWD/$out" .)
echo "wrote $out"
unzip -l "$out" | tail -1
