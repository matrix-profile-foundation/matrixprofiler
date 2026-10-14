# Computes the Distance between the 'data' and the 'query'.

This algorithm will use a rolling window, to computes the distance
thorough the whole data. This means that the minimum distance found is
the *motif* and the maximum distance is the *discord* on that time
series. **Attention** you need first to create an object using
`mass_pre()`. Read below.

## Usage

``` r
mass(
  pre_obj,
  data,
  query = data,
  index = 1,
  version = c("v3", "v2"),
  n_workers = 1
)

mass_pre(
  data,
  window_size,
  query = NULL,
  type = c("normalized", "non_normalized", "absolute", "weighted"),
  weights = NULL
)
```

## Arguments

- pre_obj:

  Required. This is the object resulting from `mass_pre()`. The is no
  *MASS* without a *pre*.

- data:

  Required. Any 1-dimension series of numbers (`matrix`, `vector`, `ts`
  etc.)

- query:

  Optional. Accepts the same types as `data` and is used for
  join-similarity. Defaults to `data` for self-similarity. **IMPORTANT**
  Use the same data used on `mass_pre()`.

- index:

  An `integer`. This is the index of the rolling window that will be
  used. Must be between `1` and `length(data) - window_size + 1`.

- version:

  A `string`. Chooses the version of MASS what will be used. Ignored if
  `mass_pre()` is not the "normalized" type.

- n_workers:

  An `integer`. The number of threads used for computing. Defaults to
  `1`. For parallel execution across processes, prefer a PSOCK cluster
  and use `n_workers = 1` within each worker. See
  [matrixprofiler](https://matrix-profile-foundation.github.io/matrixprofiler/reference/matrixprofiler.md)
  for fork limitations.

- window_size:

  Required. An integer defining the rolling window size.

- type:

  This changes how the MASS algorithm will compare the rolling window
  and the data. (See details).

- weights:

  Optional. It is used when the `type` is `weighted`, and has to be the
  same size as the `window_size`.

## Value

`mass()` returns a `list` with the `distance_profile` and the
`last_product` that is only useful for computing the Matrix Profile.

`mass_pre()` returns a `list` with several precomputations to be used on
MASS later. **Attention** use this before `mass()`.

## Details

There are currently four ways to compare the window with the data:

1.  **normalized**: this normalizes the data and the query window. This
    is the most frequently used.

2.  **non_normalized**: this won't normalize the query window. The data
    still being normalized.

3.  **absolute**: this won't normalize both the data and the query
    window.

4.  **weighted**: this normalizes the data and query window, and also
    apply a weight vector on the query.

## Examples

``` r
pre <- mass_pre(motifs_discords_small, 50)
dist_profile <- mass(pre, motifs_discords_small)
pre <- mass_pre(motifs_discords_small, 50)
dist_profile <- mass(pre, motifs_discords_small)
```
