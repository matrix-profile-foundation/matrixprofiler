# Matrix Profile Computation

STAMP Computes the best so far Matrix Profile and Profile Index for
Univariate Time Series.

STOMP is a faster implementation with the caveat that is not anytime as
STAMP or SCRIMP.

SCRIMP is a faster implementation, like STOMP, but has the ability to
return anytime results as STAMP.

MPX is by far the fastest implementation with the caveat that is not
anytime as STAMP or SCRIMP.

## Usage

``` r
stamp(
  data,
  window_size,
  query = NULL,
  exclusion_zone = 0.5,
  s_size = 1,
  n_workers = 1,
  progress = TRUE
)

stomp(
  data,
  window_size,
  query = NULL,
  exclusion_zone = 0.5,
  n_workers = 1,
  progress = TRUE,
  left_right_profile = FALSE
)

scrimp(
  data,
  window_size,
  query = NULL,
  exclusion_zone = 0.5,
  s_size = 1,
  pre_scrimp = 0.25,
  n_workers = 1,
  progress = TRUE
)

mpx(
  data,
  window_size,
  query = NULL,
  exclusion_zone = 0.5,
  s_size = 1,
  idxs = TRUE,
  distance = c("euclidean", "pearson"),
  n_workers = 1,
  progress = TRUE,
  na_strategy = c("monolithic", "segmented", "segmented_native"),
  backend = c("cpu", "cuda")
)

mpx_stream_right(
  new_data,
  start_obj,
  history = 0,
  mp_time_constraint = 0,
  progress = TRUE
)

mpx_stream_s_right(
  data,
  batch_size,
  start_obj,
  stats,
  history = 0,
  mp_time_constraint = 0,
  progress = TRUE,
  threshold = -1
)
```

## Arguments

- data:

  Required. Any 1-dimension series of numbers (`matrix`, `vector`, `ts`
  etc.) (See details).

- window_size:

  Required. An integer defining the rolling window size.

- query:

  Optional. Another 1-dimension series of numbers for an AB-join
  similarity. Default is `NULL` (See details).

- exclusion_zone:

  A numeric. Defines the size of the area around the rolling window that
  will be ignored to avoid trivial matches. Default is `0.5`, i.e., half
  of the `window_size`.

- s_size:

  A numeric. Used on anytime algorithms (stamp, scrimp, mpx) if only
  part of the computation is needed. Default is `1.0` (means 100%). For
  MPX specifically, values below `1` use a randomized subset of
  diagonals and are the recommended setting for an anytime
  approximation. At `1`, MPX computes every diagonal in monotone order
  for performance; the completed result is exact, but an interrupted
  partial result can be less representative than a randomized partial
  run.

- n_workers:

  An integer. The number of threads using for computing. Defaults to
  `1`. For parallel execution across processes, prefer a PSOCK cluster
  and use `n_workers = 1` within each worker. See
  [matrixprofiler](https://matrix-profile-foundation.github.io/matrixprofiler/reference/matrixprofiler.md)
  for fork limitations.

- progress:

  A logical. If `TRUE` (the default) will show a progress bar. Useful
  for long computations. (See details)

- left_right_profile:

  (`stomp()` only) A boolean. If `TRUE`, the function will return the
  left and right profiles.

- pre_scrimp:

  A numeric. If not zero, pre_scrimp is computed, using a fraction of
  the data. Default is `0.25`. This parameter is ignored when using
  multithread or AB-join.

- idxs:

  (`mpx()` only) A logical. Specifies if the computation will return the
  Profile Index or not. Defaults to `TRUE`.

- distance:

  (`mpx()` only) A string. Currently accepts `euclidean` and `pearson`.
  Defaults to `euclidean`.

- na_strategy:

  (`mpx()` only) Strategy for a parallel NA-aware join: `"monolithic"`
  (the default), `"segmented"`, or `"segmented_native"`. The native
  strategy supports sampled diagonal tiles when `s_size < 1`, marks the
  result as partial, and shares one native pre-calculation and scheduler
  across finite blocks. The native strategy also masks constant or
  otherwise non-normalizable windows. The `"segmented"` reference
  strategy remains restricted to inputs whose invalid windows are caused
  only by non-finite barriers.

- backend:

  (`mpx()` only) Computation backend: `"cpu"` (default) or `"cuda"`. The
  CUDA backend is optional, supports exact joins only (`s_size = 1`),
  and requires a local build with CUDA enabled.

- new_data:

  a vector with the new data

- start_obj:

  obj from another run

- history:

  an `int`. Max distance where to look for the best match in matrix
  profile. (default is NULL).

- mp_time_constraint:

  time constraint for the mp update in terms of indices. This parameter
  specifies how far back in the time series the algorithm will search
  for the nearest neighbor. A value of 0 (default) means no constraint
  and the algorithm will check the entire time series. Higher values
  limit the search to a window that extends `mp_time_constraint` indices
  into the past.

- batch_size:

  size of the new data

- stats:

  obj with stats for mpx

- threshold:

  Optional numeric threshold for distance comparison. Defaults to -1.0.

## Value

Returns a `list` with the `matrix_profile`, `profile_index` (if `idxs`
is `TRUE` in `mpx()`), and some information about the settings used to
build it, like `ez` and `partial` when the algorithm is finished early.

## Details

The Matrix Profile, has the potential to revolutionize time series data
mining because of its generality, versatility, simplicity and
scalability. In particular it has implications for time series motif
discovery, time series joins, shapelet discovery (classification),
density estimation, semantic segmentation, visualization, rule
discovery, clustering etc.

`progress`, it is really recommended to use it as feedback for long
computations. It indeed adds some (neglectable) overhead, but the
benefit of knowing that your computer is still computing is much bigger
than the seconds you may lose in the final benchmark. Native
multithreading uses RcppThread and is controlled by `n_workers`. See
[matrixprofiler](https://matrix-profile-foundation.github.io/matrixprofiler/reference/matrixprofiler.md)
for guidance on parallel execution across processes and fork
limitations. The `data` and `query` parameters will be internally
converted to a single vector using
[`as.numeric()`](https://rdrr.io/r/base/numeric.html), thus, bear in
mind that a multidimensional matrix may not work as you expect, but most
1-dimensional data types will work normally. If `query` is provided,
expect the same pre-procesment done for `data`; in addition,
`exclusion_zone` will be ignored and set to `0`. Both `data` and `query`
doesn't need to have the same size and they can be interchanged if both
are provided. The difference will be in the returning object. AB-Join
returns the Matrix Profile 'A' and 'B' i.e., the distance between a
rolling window from query to data and from data to query.

### stamp

The anytime STAMP computes the Matrix Profile and Profile Index in such
manner that it can be stopped before its complete calculation and return
the best so far results allowing ultra-fast approximate solutions.

### stomp

The STOMP uses a faster implementation to compute the Matrix Profile and
Profile Index. It can be stopped earlier by the user, but the result is
not considered anytime, just incomplete. For a anytime algorithm, use
`stamp()` or `scrimp()`.

### scrimp

The SCRIMP algorithm was the anytime solution for stomp. It is as fast
as stomp but allows the user to cancel the computation and get an
approximation of the final result. The exact SCRIMP phase shares the
NA-aware diagonal recurrence used by MPX. In serial mode, pre-scrimp
first seeds the profile from sampled distance profiles; the
multithreaded implementation skips this pre-scrimp stage.

### mpx

This algorithm was developed apart from the main Matrix Profile branch
that relies on Fast Fourier Transform (FFT) at least in one part of the
process. This algorithm doesn't use FFT at all and is several times
faster. It also relies on Ogita's work for better precision computing
mean and standard deviation (part of the process).

Inputs containing `NA`, `NaN`, `Inf`, or `-Inf` are automatically routed
to the missing-data implementation. Fully finite inputs continue to use
the faster regular MPX implementation.

For MPX, values of `s_size` below `1` select a randomized subset of
diagonals and provide the intended anytime approximation. With
`s_size = 1`, all diagonals are processed in monotone order for better
locality, so the final result is exact but an interruption before
completion can be less representative than a randomized partial run. Use
`s_size < 1` when an interrupted result should remain an anytime
approximation.

The native segmented parallel AA implementation has an optional
cache-tuning setting for large NA-aware self-joins. This setting applies
only when `na_strategy = "segmented_native"` selects the native NA-aware
parallel path. When `s_size = 1`, `MATRIXPROFILER_NATIVE_AA_BLOCK_SIZE`
controls the number of profile windows in each rectangular internal
block. Its default is `8192`; values from `1024` to `262144` are
accepted and invalid values fall back to the default. For example, set
it for one shell invocation with
`MATRIXPROFILER_NATIVE_AA_BLOCK_SIZE=16384 Rscript script.R`, or set it
in R with `Sys.setenv(MATRIXPROFILER_NATIVE_AA_BLOCK_SIZE = "16384")`
before calling `mpx()`. Remove the override with
`Sys.unsetenv("MATRIXPROFILER_NATIVE_AA_BLOCK_SIZE")`. The setting is
read at call time, does not change the public `mpx()` API or its
numerical parameters, and is ignored for AB joins, finite inputs using
regular MPX, and sampled native runs with `s_size < 1`. It primarily
trades cache locality against task overhead; `8192` is a suitable
starting point, while larger values should only be selected after
benchmarking a substantially larger workload.

The internal C++ target `matrix_profile_native_pairs_per_task`
(currently `100000000` pairwise window comparisons) controls the
granularity of native parallel AB tasks and sampled AA tasks. It is used
only by the native NA-aware parallel implementation. It is not an
`mpx()` argument or an environment variable; changing it requires a
source change and recompiling the package. It is normally left
unchanged: increase it only when profiling a very large workload shows
excessive scheduling overhead, while retaining several tasks per worker
for load balancing. Neither setting affects the regular MPX paths or the
serial NA-aware implementations.

## This document

Last updated on 2026-10-14 using R version 4.6.1.

## References

- Yeh CCM, Zhu Y, Ulanova L, Begum N, Ding Y, Dau HA, et al. Matrix
  profile I: All pairs similarity joins for time series: A unifying view
  that includes motifs, discords and shapelets. Proc - IEEE Int Conf
  Data Mining, ICDM. 2017;1317-22.

&nbsp;

- Zhu Y, Imamura m, Nikovski D, Keogh E. Matrix Profile VII: Time Series
  Chains: A New Primitive for Time Series Data Mining. Knowl Inf Syst.
  2018 Jun 2;1-27.

&nbsp;

- Zhu Y, Zimmerman Z, Senobari NS, Yeh CM, Funning G. Matrix Profile II
  : Exploiting a Novel Algorithm and GPUs to Break the One Hundred
  Million Barrier for Time Series Motifs and Joins. Icdm. 2016 Jan
  22;54(1):739-48.

Website:
[http://www.cs.ucr.edu/~eamonn/MatrixProfile.html](http://www.cs.ucr.edu/~eamonn/MatrixProfile.md)

## See also

[`mass()`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/mass.md)
for the underlying algorithm that finds best match of a query.

`mpxab()` for the forward and reverse join-similarity.

Other matrix profile computations:
[`mstomp()`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/mstomp.md)

## Examples

``` r
mp <- stamp(motifs_discords_small, 50)
#> Computing: [========================================] 100% (done)                         
mp <- stomp(motifs_discords_small, 50)
mp <- scrimp(motifs_discords_small, 50)
#> Computing: [========================================] 100% (done)                         
#> Computing: [========================================] 100% (done)                         
mp <- mpx(motifs_discords_small, 50)
#> Computing: [========================================] 100% (done)                         
```
