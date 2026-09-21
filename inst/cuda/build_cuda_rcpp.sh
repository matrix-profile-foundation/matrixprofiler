#!/usr/bin/env sh
# Build matrixprofiler locally with its optional CUDA MPX proof of concept.
# This script never changes src/Makevars and must not be used for CRAN builds.

set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cuda_root=$(CDPATH= cd -- "$(dirname -- "$(command -v nvcc)")/.." && pwd)
object="${TMPDIR:-/tmp}/matrixprofiler_mpx_gpu_poc.o"

nvcc -std=c++17 -O3 -arch=sm_89 -Xcompiler -fPIC -c "$root/inst/cuda/mpx_gpu_poc.cu" -o "$object"

export MATRIXPROFILER_ENABLE_CUDA=1
export MATRIXPROFILER_CUDA_OBJECT="$object"
export MATRIXPROFILER_CUDA_LIBS="-L$cuda_root/lib64 -lcudart"

cd "$root"
Rscript -e 'renv::load(); pkgbuild::clean_dll(); devtools::load_all(debug = FALSE, quiet = FALSE)'

if [ "$#" -gt 0 ]; then
	exec "$@"
fi