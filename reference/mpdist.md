# MPdist - Distance between Time Series using Matrix Profile

MPdist is a distance measure which considers two time series to be
similar if they share many similar subsequences, regardless of the order
of matching subsequences.

## Usage

``` r
mpdist(
  ref_data,
  query_data,
  window_size,
  type = c("simple", "vector"),
  thr = 0.05,
  n_workers = 1L
)
```

## Arguments

- ref_data:

  a `matrix` or a `vector`. The reference data.

- query_data:

  a `matrix` or a `vector`. The query data.

- window_size:

  an `int`. Size of the sliding window. Must be at least 4.

- type:

  the type of result: `"simple"` (default) returns a scalar distance,
  `"vector"` returns distance between all sliding windows.

- thr:

  threshold for MPdist. (Default is `0.05`).

- n_workers:

  an `int`. Number of threads for parallel computation. (Default is
  `1`). For parallel execution across processes, prefer a PSOCK cluster
  and use `n_workers = 1` within each worker. See
  [matrixprofiler](https://matrix-profile-foundation.github.io/matrixprofiler/reference/matrixprofiler.md)
  for fork limitations.

## Value

Returns the distance scalar (`"simple"`) or a numeric vector
(`"vector"`).

## References

- Gharghabi S, Imani S, Bagnall A, Darvishzadeh A, Keogh E. Matrix
  Profile XII: MPdist: A Novel Time Series Distance Measure to Allow
  Data Mining in More Challenging Scenarios. In: 2018 IEEE International
  Conference on Data Mining (ICDM). 2018.

Website: <https://sites.google.com/site/mpdistinfo/>

## Examples

``` r
ref_data <- motifs_discords_small[1:200]
qe_data <- motifs_discords_small[201:400]
w <- 30
deq <- mpdist(ref_data, qe_data, w)
```
