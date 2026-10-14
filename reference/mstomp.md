# Multivariate STOMP Algorithm

Computes the Matrix Profile and Profile Index for Multivariate Time
Series.

## Usage

``` r
mstomp(
  data,
  window_size,
  exclusion_zone = 0.5,
  must_dim = NULL,
  exc_dim = NULL,
  n_workers = 1L,
  progress = TRUE
)
```

## Arguments

- data:

  a `matrix` of `numeric`, where each column is a time series. Accepts
  `vector`, `list` and `data.frame` too.

- window_size:

  an `int` with the size of the sliding window.

- exclusion_zone:

  a `numeric`. Size of the exclusion zone, based on window size (default
  is `0.5`).

- must_dim:

  an `int` or `vector` of which dimensions to forcibly include (default
  is `NULL`).

- exc_dim:

  an `int` or `vector` of which dimensions to exclude (default is
  `NULL`).

- n_workers:

  an `int`. Number of workers for parallel execution (default is `1`).

- progress:

  a `logical`. If `TRUE` will show a progress bar. (Default is `TRUE`).

## Value

Returns a `list` with the matrix profile `matrix_profile`, profile index
`profile_index`, left and right matrix profiles `left_matrix_profile`,
`right_matrix_profile` and profile indices `left_profile_index`,
`right_profile_index`, window size `w`, number of dimensions `n_dim`,
exclusion zone `ez`, must dimensions `must` and excluded dimensions
`exc`.

## References

- Yeh CM, Kavantzas N, Keogh E. Matrix Profile VI : Meaningful
  Multidimensional Motif Discovery.

&nbsp;

- Zhu Y, Imamura M, Nikovski D, Keogh E. Matrix Profile VII: Time Series
  Chains: A New Primitive for Time Series Data Mining. Knowl Inf Syst.
  2018 Jun 2;1-27.

Website: <https://sites.google.com/view/mstamp/>

Website:
[http://www.cs.ucr.edu/~eamonn/MatrixProfile.html](http://www.cs.ucr.edu/~eamonn/MatrixProfile.md)

## See also

Other matrix profile computations:
[`stamp()`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/mp_algos.md)

## Examples

``` r
mp <- mstomp(motifs_discords_small, 50, progress = FALSE)
```
