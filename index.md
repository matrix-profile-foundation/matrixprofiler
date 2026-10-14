# README

Francisco Bischoff

14 Oct 2026

![matrixprofiler logo](reference/figures/logo.png)

## Matrix Profile for R

## What’s New — v0.1.12 (October 2026)

> **`matrixprofiler` is the high-performance C++ engine** of the Matrix
> Profile ecosystem. This release focuses on native parallelism, new
> algorithms, and safer NA handling.

Key highlights:

- ⚡ **`RcppThread`** replaces `RcppParallel`/Intel TBB —
  zero-dependency, fork-safe multithreading on all platforms.
- 🆕
  **[`contrast()`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/contrast.md)**
  — native Contrast Profile computation with multi-threading, Pearson
  and Euclidean distance metrics.
- 🆕
  **[`rfcp()`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/rfcp.md)
  &
  [`rfcp_top_k()`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/rfcp_top_k.md)**
  — Relative Frequency Contrast Profile with Top-K morphology discovery,
  frequency adaptation, checkpoints, and CUDA acceleration.
- 🆕
  **[`mstomp()`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/mstomp.md)**
  — native C++ Multivariate STOMP with `n_workers` support.
- 🆕
  **[`mpdist()`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/mpdist.md)**
  — native C++ MPdist (scalar and vector variants), multi-thread.
- 🔢 **`NA_integer_`** — invalid profile indices now follow R’s
  canonical missing-value convention (was `-1L`).
- 🛡️ **`na_strategy`** for all `mov_*` window functions — predictable,
  non-cascading NA handling.

📄 Full API changes, breaking changes, and migration guide:
**[MIGRATION.md](https://github.com/matrix-profile-foundation/matrixprofiler/blob/main/MIGRATION.md)**

### Overview

R package implementing the UCR Matrix Profile Algorithm
([http://www.cs.ucr.edu/~eamonn/MatrixProfile.html](http://www.cs.ucr.edu/~eamonn/MatrixProfile.md))
as an optimised C++17 engine.

This package provides the **core computational primitives** for building
and querying Matrix Profiles. Higher-level analytics (motif search,
segmentation, snippets, plotting) live in the
[`tsmp`](https://CRAN.R-project.org/package=tsmp) package, which uses
`matrixprofiler` as its backend.

### Installation

``` r

# Install the released version from CRAN
install.packages("matrixprofiler")
# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("matrix-profile-foundation/matrixprofiler")
```

#### Optional: NVIDIA CUDA GPU Acceleration

By default, `matrixprofiler` installs as a portable, high-performance
CPU engine with multi-threading via `RcppThread` (zero external
dependencies, 100% CRAN compliant).

If you have an **NVIDIA GPU** and the **CUDA Toolkit** installed, you
can compile `matrixprofiler` with native CUDA acceleration for supported
algorithms (such as
[`rfcp()`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/rfcp.md)
and
[`mpx()`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/mp_algos.md)):

``` bash
# Clone the repository
git clone https://github.com/matrix-profile-foundation/matrixprofiler.git
cd matrixprofiler

# Compile CUDA kernels and install directly into your R library
./inst/cuda/build_cuda_rcpp.sh R CMD INSTALL .
# or with devtools:
./inst/cuda/build_cuda_rcpp.sh Rscript -e 'devtools::install()'
```

Once installed, check if CUDA is enabled in your R session:

``` r

library(matrixprofiler)
has_cuda()
#> [1] TRUE
```

### Currently available Features

- **STAMP** (single and multi-thread)
- **STOMP** (single and multi-thread)
- **SCRIMP** (single and multi-thread, not for AB-joins yet)
- **MPX** (single and multi-thread)
- **Contrast Profile (`contrast`)** — Contrast profile between positive
  and negative classes (Euclidean/Pearson, multi-thread)
- **RFCP (`rfcp`)** — Relative Frequency Contrast Profile for
  discriminative pattern discovery (CPU multi-threaded & CUDA GPU)
- **Top-K Morphology Discovery (`rfcp_top_k`)** — Iterative discovery of
  top-K distinct morphologies with frequency adaptation and
  checkpointing
- **mSTOMP** — Multivariate STOMP, native C++ with `n_workers`
- **MPdist** — Matrix Profile Distance, scalar and vector variants,
  multi-thread
- Misc:
  - MASS v2.0 / v3.0 (multi-variant)
  - MASS extensions: WQ (Weighted Query), ABS (Absolute), UN
    (Unnormalized)
  - Sliding-window functions (`mov_mean`, `mov_std`, `mov_var`,
    `mov_min`, `mov_max`, `muinvn`, `zero_crossing`) with
    `na_strategy = c("window_na", "remove", "error")` and `na.rm`
    control

### Roadmap

Replace Intel TBB / `RcppParallel` with `RcppThread` (fork-safe, zero
external dep)

Native Multivariate STOMP
([`mstomp()`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/mstomp.md))
in C++

Native MPdist
([`mpdist()`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/mpdist.md))
in C++

Contrast Profile
([`contrast()`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/contrast.md))
in native C++

Relative Frequency Contrast Profile
([`rfcp()`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/rfcp.md))
and Top-K Morphology Discovery
([`rfcp_top_k()`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/rfcp_top_k.md))

Standardise invalid indices to `NA_integer_`

`na_strategy` for all sliding-window functions

STOMPi (On-line / streaming version)

SiMPle-Fast (Music Analysis)

VALMOD (Variable Length Motifs) — candidate for `matrixprofiler` core

MASS extensions: ADP (Approximate Distance Profile, with PAA), QwG
(Query with Gap)

### Other projects with Matrix Profile

- Python: <https://github.com/target/matrixprofile-ts>
- Python: <https://github.com/ZiyaoWei/pyMatrixProfile>
- Python: <https://github.com/juanbeleno/owlpy>
- Python: <https://github.com/javidlakha/matrix-profile>
- Python: <https://github.com/shapelets/khiva-python>
- R: <https://github.com/shapelets/khiva-r>
- Matlab: <https://github.com/shapelets/khiva-matlab>
- Java: <https://github.com/shapelets/khiva-java>
- Java: <https://github.com/ensozos/Matrix-Profile>
- Kotlin: <https://github.com/shapelets/khiva-kotlin>
- C++ (CUDA and OPENCL): <https://github.com/shapelets/khiva>
- CUDA: <https://github.com/zpzim/STOMPSelfJoin>
- CUDA: <https://github.com/zpzim/SCAMP>

### Matrix Profile Foundation

Our next step unifying the Matrix Profile implementation in several
programming languages.

Visit: [Matrix Profile
Foundation](https://matrix-profile-foundation.github.io/)

### Benchmarks

Available at [RPubs](https://rpubs.com/franzbischoff/matrixprofiler).

### Package dependencies

![Matrixprofiler package dependencies network
graph](reference/figures/dependency_plot-1.png)

### Code of Conduct

Please note that the matrixprofiler project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
