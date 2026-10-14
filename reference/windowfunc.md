# Several moving window functions

These functions compute running window statistics with comprehensive
support for missing and non-finite values (`NA`, `NaN`, `Inf`).

## Usage

``` r
mov_mean(
  data,
  window_size,
  type = c("ogita", "normal", "weighted", "fading"),
  eps = 0.9,
  na_strategy = c("window_na", "remove", "error"),
  na.rm = NULL
)

mov_var(
  data,
  window_size,
  type = c("ogita", "normal", "weighted", "fading"),
  eps = 0.9,
  na_strategy = c("window_na", "remove", "error"),
  na.rm = NULL
)

mov_sum(
  data,
  window_size,
  type = c("ogita", "normal", "weighted", "fading"),
  eps = 0.9,
  na_strategy = c("window_na", "remove", "error"),
  na.rm = NULL
)

mov_max(
  data,
  window_size,
  na_strategy = c("window_na", "remove", "error"),
  na.rm = NULL
)

mov_min(
  data,
  window_size,
  na_strategy = c("window_na", "remove", "error"),
  na.rm = NULL
)

mov_std(
  data,
  window_size,
  rcpp = TRUE,
  na_strategy = c("window_na", "remove", "error"),
  na.rm = NULL
)

movmean_std(
  data,
  window_size,
  rcpp = TRUE,
  na_strategy = c("window_na", "remove", "error"),
  na.rm = NULL
)

muinvn(
  data,
  window_size,
  n_workers = 1L,
  na_strategy = c("window_na", "remove", "error"),
  na.rm = NULL
)

zero_crossing(data, window_size)
```

## Arguments

- data:

  A `vector` or a column `matrix` of `numeric`.

- window_size:

  An `integer`. The size of the rolling window.

- type:

  A `string`. Select between several algorithms. Default is `ogita` (See
  details).

- eps:

  A `numeric`. Used only for fading algorithms (See details), otherwise
  has no effect.

- na_strategy:

  A `string` indicating how missing values (`NA`/`NaN`) are handled:

  - `"window_na"` (default): Any window touching an `NA`/`NaN` yields
    `NA_real_`, without contaminating adjacent complete windows or
    stopping with an error.

  - `"remove"`: Missing values inside the window are ignored, computing
    the statistic over remaining finite values.

  - `"error"`: Throws an error immediately if missing values are
    present.

- na.rm:

  A `logical`. Standard R compatibility shortcut. If provided, `TRUE`
  maps to `na_strategy = "remove"`, and `FALSE` maps to
  `na_strategy = "window_na"`.

- rcpp:

  A `logical`. If `TRUE` will use the Rcpp implementation, otherwise
  will use the R implementation, that may or not be slower.

- n_workers:

  An `integer`. The number of threads using for computing. Defaults to
  `1`. For parallel execution across processes, prefer a PSOCK cluster
  and use `n_workers = 1` within each worker. See
  [matrixprofiler](https://matrix-profile-foundation.github.io/matrixprofiler/reference/matrixprofiler.md)
  for fork limitations.

## Value

`mov_mean()` returns a `vector` with moving `avg`.

`mov_var()` returns a `vector` with moving `var`.

`mov_sum()` returns a `vector` with moving `sum`.

`mov_max()` returns a `vector` with moving `max`.

`mov_min()` returns a `vector` with moving `min`.

`mov_std()` returns a `vector` with moving `sd`.

`movmean_std()` returns a `list` with `vectors` of the moving `avg`,
`sd`, `sig`, `sum` and `sqrsum`.

`muinvn()` returns a `list` with `vectors` of moving `avg` and `sig`.

`zero_crossing()` returns a `vector` of times the data crossed the
'zero' line inside a rolling window.

## Details

Some functions may use different algorithms to compute the results. The
available types are:

1.  **ogita**: This is the default. It uses the Ogita *et al.*, Accurate
    Sum, and Dot Product for precision. It is not the fastest algorithm,
    but the time spent vs. guarantee of precision worth it.

2.  **normal**: This uses the `cumsum` method that is faster, but
    unreliable in some situations (I have to find the references, but is
    true).

3.  **weighted**: This uses Rodrigues P., *et al.* algorithm that uses a
    weighted window for online purposes. The `eps` argument controls the
    factor. (The function is not online yet)

4.  **fading**: This also uses Rodrigues P., *et al.* algorithm that in
    this case, uses a fading factor, also for online purposes. he `eps`
    argument controls the factor. (The function is not online yet)

Another important detail is that the *standard deviation* we use for all
computations is the *population* (i.e.: divided by `n`), not the
*sample* (i.e.: divided by `n - 1`). That is why we also provide the
internally the `:::std()` function that computes the *population*,
differently from [`stats::sd()`](https://rdrr.io/r/stats/sd.html) that
is the *sample* kind. Further more, `movmean_std()` shall be used when
you need both results in one computation. This is faster than call
`mov_mean()` followed by `mov_std()`. Finally, `muinvn()` is kept like
that for historical reasons, as it is the function used by
[`mpx()`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/mp_algos.md).
It returns the `sig` (stable inverse centered norm) instead of `std`
(`sig` is equals to `1 / (std * sqrt(window_size))`).

## Examples

``` r
mov <- mov_mean(motifs_discords_small, 50)
mov <- mov_var(motifs_discords_small, 50)
mov <- mov_sum(motifs_discords_small, 50)
mov <- mov_max(motifs_discords_small, 50)
mov <- mov_min(motifs_discords_small, 50)
mov <- mov_std(motifs_discords_small, 50)
mov <- movmean_std(motifs_discords_small, 50)
mov <- muinvn(motifs_discords_small, 50)
zero_cross <- zero_crossing(motifs_discords_small, 50)
```
