# Relative Frequency Contrast Profile

Computes an exact, NA-aware Relative Frequency Contrast Profile (RFCP)
between a positive series and a negative series. For every valid
positive window, RFCP compares ordered non-overlapping neighbor ranks
from the positive self-join and the positive-to-negative AB-join.

## Usage

``` r
rfcp(
  positive_data,
  negative_data,
  window_size,
  max_freq,
  exclusion_radius = NULL,
  query_begin = 1L,
  query_end = NULL,
  return_profiles = FALSE,
  n_workers = 1L,
  progress = TRUE,
  backend = "cpu"
)
```

## Arguments

- positive_data:

  Numeric positive/query time series.

- negative_data:

  Numeric negative/reference time series.

- window_size:

  Integer subsequence length, at least `2` and no larger than either
  input series.

- max_freq:

  Integer number of non-overlapping neighbor ranks requested.

- exclusion_radius:

  Integer exclusion radius between accepted neighbor starts. `NULL` uses
  `window_size`, the RFCP protocol setting.

- query_begin:

  First positive query start to process, using one-based inclusive
  coordinates. Defaults to the first start.

- query_end:

  Last positive query start to process, using one-based inclusive
  coordinates. `NULL` processes through the final start.

- return_profiles:

  Whether to return the complete RFMP, RFCP, and index matrices for the
  requested query block. This validation mode is limited to 5,000,000
  rank-query cells. Defaults to `FALSE`.

- n_workers:

  Number of native worker threads for the CPU backend. Ignored when
  `backend = "cuda"`. Defaults to `1`.

- progress:

  Logical progress flag. Reserved for native progress reporting;
  defaults to `TRUE`.

- backend:

  Computation backend: `"cpu"` (default) or `"cuda"`. The CUDA backend
  is optional and requires a local build with CUDA enabled (see
  `inst/cuda/build_cuda_rcpp.sh`).

## Value

A list with RFCP outputs:

- `rms_profile`: numeric vector of root-mean-square contrasts.

- `valid_window_positive`: logical mask of valid query windows.

- `valid_window_negative`: logical mask of valid reference windows.

- `plato_index`: one-based index of the highest-contrast query.

- `plato_rms`: contrast score at `plato_index`.

- `plato_twin_index`: one-based index of the nearest self-join neighbor
  at `plato_index`.

- `rfcp_at_plato`: clipped contrasts across accepted neighbor ranks.

- `rfmp_aa_at_plato`: positive self-join distances at plato.

- `rfmp_ab_at_plato`: negative AB distances at plato.

- `rfmpi_aa_at_plato`: positive neighbor indices at plato.

- `rfmpi_ab_at_plato`: negative neighbor indices at plato.

- `requested_max_freq`: requested rank count.

- `effective_max_freq`: highest paired rank observed.

- `query_begin`, `query_end`: one-based query range.

- `completed_query_end`: end of completed queries.

- `partial`: whether computation stopped early.

- `rfmp_aa`, `rfmp_ab`, `rfcp`, `rfmpi_aa`, `rfmpi_ab`: optional full
  profile matrices when `return_profiles = TRUE`.

## Details

`NA`, `NaN`, `Inf`, and `-Inf` are barriers. Windows touching a barrier,
constant windows, and otherwise non-normalizable windows cannot be
queries or neighbors. Returned neighbor indices are one-based global
coordinates. Candidates are ordered by increasing distance and then
increasing index. After accepting a neighbor at `j`, the closed interval
`[j - exclusion_radius, j + exclusion_radius]` is excluded. The same
neighborhood around the query is excluded before AA selection.

For rank `k`, distances are clipped at `sqrt(2 * window_size)` before
calculating `max((AB - AA) / sqrt(2 * window_size), 0)`. Missing ranks
contribute zero while the RMSC denominator remains the requested
`max_freq`. `effective_max_freq` is diagnostic: it is the largest number
of paired finite AA/AB ranks observed for any completed query and never
changes the denominator.

The default mode materializes only the RMSC values for the requested
query block and the rank vectors for the selected Plato. Blocks use
global query coordinates and can be concatenated for checkpointed
execution. With `return_profiles = TRUE`, RFMP matrices contain
unclipped z-normalized Euclidean distances and RFCP contains the
corresponding clipped contrasts. Native query tiles use
`min(32768, max(1024, 16 * window_size))`, further capped at the ceiling
of the requested query count divided by the effective worker count. The
first tile starts at `query_begin`, so checkpoint blocks do not perform
discarded warm-up rows. The CUDA backend executes the full join
query-by-query with optional profile matrices returned identically to
the CPU backend.

## Examples

``` r
pos <- c(sin(seq(0, 4 * pi, length.out = 100)))
neg <- c(cos(seq(0, 10 * pi, length.out = 150)))
res <- rfcp(pos, neg, window_size = 10, max_freq = 2)
res$plato_index
#> [1] 20
```
