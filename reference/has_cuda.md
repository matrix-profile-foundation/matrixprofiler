# Check if package was compiled with CUDA support

Queries whether `matrixprofiler` was compiled with NVIDIA CUDA GPU
acceleration (`MATRIXPROFILER_ENABLE_CUDA=1`).

## Usage

``` r
has_cuda()
```

## Value

A logical value: `TRUE` if CUDA support was enabled at compile time,
`FALSE` otherwise.

## Details

By default, `matrixprofiler` is compiled without CUDA to maintain
portability and compliance with CRAN policies, running all algorithms on
multi-threaded CPU engines powered by `RcppThread`.

If an NVIDIA GPU and the CUDA Toolkit are available on the host system,
the package can be recompiled with CUDA acceleration for supported
algorithms (such as
[`rfcp`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/rfcp.md)
and
[`mpx`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/mp_algos.md))
using the helper script `inst/cuda/build_cuda_rcpp.sh` or by setting the
environment variables:

- `MATRIXPROFILER_ENABLE_CUDA=1`

- `MATRIXPROFILER_CUDA_OBJECT` pointing to precompiled CUDA object files
  (`.o`)

- `MATRIXPROFILER_CUDA_LIBS` pointing to the CUDA runtime library (e.g.
  `-lcudart`)

When `has_cuda()` returns `FALSE`, attempting to execute an algorithm
with `backend = "cuda"` will stop with an informative error instructing
how to enable CUDA.

## Examples

``` r
has_cuda()
#> [1] FALSE
```
