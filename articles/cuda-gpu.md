# NVIDIA CUDA GPU Acceleration

## Overview

By default, **`matrixprofiler`** compiles as a portable,
high-performance CPU package using **`RcppThread`** for multi-threading
across Linux, macOS, and Windows. This ensures 100% compliance with CRAN
policies without requiring proprietary drivers or external hardware
dependencies.

For users with **NVIDIA GPUs**, `matrixprofiler` provides native CUDA
backends for computationally demanding tasks:

- **`rfcp(..., backend = "cuda")`**: Relative Frequency Contrast
  Profile. Accelerated candidate sorting, greedy exclusion zone
  suppression, and Top-K Plato discovery completely on the GPU with up
  to 55x speedup over CPU-sorting GPU baselines and over 12x speedup
  compared to single-threaded CPU execution.
- **`mpx(..., backend = "cuda")`**: GPU-accelerated Matrix Profile for
  continuous time series.

------------------------------------------------------------------------

## Hardware and Software Prerequisites

To build and run `matrixprofiler` with CUDA acceleration, you need:

1.  **NVIDIA GPU** with Compute Capability \>= 6.0 (Pascal, Volta,
    Turing, Ampere, Ada Lovelace, Hopper, Blackwell).

2.  **NVIDIA CUDA Toolkit** (version 11.0 or newer) with `nvcc`
    accessible in your system `$PATH`:

    ``` bash
    nvcc --version
    nvidia-smi
    ```

3.  A standard C++17 compiler (`g++` or `clang++`).

------------------------------------------------------------------------

## Building matrixprofiler with CUDA

### Method 1: Using the Included Build Script (Recommended)

From the cloned repository (or unpacked source package), run the
provided helper script
[`inst/cuda/build_cuda_rcpp.sh`](https://github.com/matrix-profile-foundation/matrixprofiler/blob/main/inst/cuda/build_cuda_rcpp.sh):

``` bash
# Clone repository
git clone https://github.com/matrix-profile-foundation/matrixprofiler.git
cd matrixprofiler

# Build CUDA objects and install into your R library:
./inst/cuda/build_cuda_rcpp.sh R CMD INSTALL .

# Or if using devtools:
./inst/cuda/build_cuda_rcpp.sh Rscript -e 'devtools::install()'
```

### Method 2: Manual Environment Configuration

If you prefer to compile manually or integrate with custom build
pipelines:

``` bash
# 1. Compile CUDA objects (.cu -> .o) with position-independent code (-fPIC)
CUDA_HOME=${CUDA_HOME:-/usr/local/cuda}
nvcc -std=c++17 -O3 -Xcompiler -fPIC \
     -c inst/cuda/mpx_gpu_poc.cu -o /tmp/matrixprofiler_mpx_gpu_poc.o
nvcc -std=c++17 -O3 -Xcompiler -fPIC \
     -c inst/cuda/rfcp_gpu_poc.cu -o /tmp/matrixprofiler_rfcp_gpu_poc.o

# 2. Export build flags
export MATRIXPROFILER_ENABLE_CUDA=1
export MATRIXPROFILER_CUDA_OBJECT="/tmp/matrixprofiler_mpx_gpu_poc.o /tmp/matrixprofiler_rfcp_gpu_poc.o"
export MATRIXPROFILER_CUDA_LIBS="-L$CUDA_HOME/lib64 -lcudart"

# 3. Install the package
R CMD INSTALL .
```

------------------------------------------------------------------------

## Verifying CUDA Acceleration in R

In your R session, call
[`has_cuda()`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/has_cuda.md)
to verify that the active binary was compiled with CUDA support:

``` r

library(matrixprofiler)

has_cuda()
#> [1] TRUE
```

If
[`has_cuda()`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/has_cuda.md)
returns `FALSE`, the package is running on the standard CPU backend.

------------------------------------------------------------------------

## Usage Examples

### GPU-Accelerated RFCP (Relative Frequency Contrast Profile)

``` r

library(matrixprofiler)

# Generate synthetic series
pos_series <- cumsum(rnorm(5000))
neg_series <- cumsum(rnorm(5000))

# Execute on GPU:
if (has_cuda()) {
  result <- rfcp(
    pos_data = pos_series,
    neg_data = neg_series,
    window_size = 100,
    backend = "cuda"
  )
  print(result)
} else {
  message("CUDA is not enabled on this installation. Running on multi-core CPU...")
  result <- rfcp(pos_series, neg_series, window_size = 100, n_workers = 4)
}
```

### Top-K Plato Discovery with GPU Acceleration

``` r

if (has_cuda()) {
  top_k <- rfcp_top_k(
    pos_data = pos_series,
    neg_data = neg_series,
    window_size = 100,
    k = 5,
    backend = "cuda"
  )
  print(top_k)
}
```
