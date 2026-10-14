# Contrast Profile Computation

Computes the contrast profile of two (classes of) time series.

## Usage

``` r
contrast(
  negative_data,
  positive_data,
  window_size,
  positive_matrix = NULL,
  exclusion_zone = 0.5,
  distance = c("euclidean", "pearson"),
  n_workers = 1L,
  progress = TRUE
)
```

## Arguments

- negative_data:

  Required. Any 1-dimension series of numbers (`matrix`, `vector`, `ts`
  etc.) where the pattern is not present

- positive_data:

  Required. Any 1-dimension series of numbers (`matrix`, `vector`, `ts`
  etc.) where the pattern is present

- window_size:

  Required. An integer defining the rolling window size.

- positive_matrix:

  Optional. A precomputed self-similar matrix profile of the positive
  data.

- exclusion_zone:

  A numeric. Defines the size of the area around the rolling window that
  will be ignored to avoid trivial matches. Default is `0.5`, i.e., half
  of the `window_size`.

- distance:

  A string. Currently accepts `euclidean` and `pearson`. Defaults to
  `euclidean`.

- n_workers:

  An integer. The number of threads using for computing. Defaults to
  `1`.

- progress:

  A logical. If `TRUE` (the default) will show a progress bar. Useful
  for long computations.

## Value

Returns a `list` with the `contrast_profile`, `plato`, `plato_nn`,
`plato_idx`, `plato_nn_idx`, `w`, `ez`, `euclidean` values

## Details

### Constrast Profile

This algorithm returns the contrast profile of two time series, which
shows the position of patters that are similar in the positive data, but
at the same time very dissimilar in the negative data. In other words,
this means that such a pattern represents well positive data and may be
taken as a "signature" of that class. More information can be found in
the references.

## References

- R. Mercer, S. Alaee, A. Abdoli, S. Singh, A. Murillo and E. Keogh,
  "Matrix Profile XXIII: Contrast Profile: A Novel Time Series Primitive
  that Allows Real World Classification," 2021 IEEE International
  Conference on Data Mining (ICDM), 2021, pp. 1240-1245, doi:
  10.1109/ICDM51629.2021.00151.

Website:
[http://www.cs.ucr.edu/~eamonn/MatrixProfile.html](http://www.cs.ucr.edu/~eamonn/MatrixProfile.md)

## Examples

``` r
cp <- contrast(motifs_discords_small, rev(motifs_discords_small), 50)
#> Computing: [========================================] 100% (done)                         
```
