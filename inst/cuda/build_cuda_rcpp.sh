#!/usr/bin/env sh
# Compile matrixprofiler's optional CUDA objects and run a caller-provided
# command with the required build environment. This script never changes
# src/Makevars and must not be used for CRAN builds.

set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cuda_root=$(CDPATH= cd -- "$(dirname -- "$(command -v nvcc)")/.." && pwd)
object_mpx="${TMPDIR:-/tmp}/matrixprofiler_mpx_gpu_poc.o"
object_rfcp="${TMPDIR:-/tmp}/matrixprofiler_rfcp_gpu_poc.o"

nvcc -std=c++17 -O3 -arch=sm_89 -Xcompiler -fPIC \
     -c "$root/inst/cuda/mpx_gpu_poc.cu"  -o "$object_mpx"
nvcc -std=c++17 -O3 -arch=sm_89 -Xcompiler -fPIC \
     -c "$root/inst/cuda/rfcp_gpu_poc.cu" -o "$object_rfcp"

export MATRIXPROFILER_ENABLE_CUDA=1
export MATRIXPROFILER_CUDA_OBJECT="$object_mpx $object_rfcp"
export MATRIXPROFILER_CUDA_LIBS="-L$cuda_root/lib64 -lcudart"

if [ "$#" -eq 0 ]; then
	printf '%s\n' "Usage: $0 COMMAND [ARGUMENTS...]" >&2
	printf '%s\n' "Example: $0 Rscript -e 'renv::install(\"$root\")'" >&2
	exit 2
fi

exec "$@"