# Math Functions

`znorm()`: Normalizes data for mean Zero and Standard Deviation One

`ed_corr()`: Converts euclidean distances into correlation values

`corr_ed()`: Converts correlation values into euclidean distances

`mode()`: Returns the most common value from a vector of integers

`std()`: Population SD, as R always calculate with n-1 (sample), here we
fix it.

`normalize()`: Normalizes data to be between min and max.

`complexity()`: Computes the complexity index of the data

`binary_split()`: Creates a vector with the indexes of binary split.

## Usage

``` r
znorm(data, rcpp = TRUE)

ed_corr(data, w, rcpp = TRUE)

corr_ed(data, w, rcpp = TRUE)

mode(x, rcpp = FALSE)

std(data, na.rm = FALSE, rcpp = TRUE)

normalize(data, min_lim = 0, max_lim = 1, rcpp = FALSE)

complexity(data)

binary_split(n, rcpp = TRUE)
```

## Arguments

- data:

  a `vector` of `numeric`.

- rcpp:

  A `logical`. If `TRUE` will use the Rcpp implementation, otherwise
  will use the R implementation, that may or not be slower.

- w:

  the window size

- x:

  a `vector` of `integers`.

- na.rm:

  A logical. If `TRUE` remove the `NA` values from the computation.

- min_lim:

  A number

- max_lim:

  A number

- n:

  size of the vector

## Value

`znorm()`: Returns the normalized data

`ed_corr()`: Returns the converted values from euclidean distance to
correlation values.

`corr_ed()`: Returns the converted values from euclidean distance to
correlation values.

`mode()`: Returns the most common value from a vector of integers.

`std()`: Returns the corrected standard deviation from sample to
population.

`normalize()`: Returns the normalized data between min and max.

`complexity()`: Returns the complexity index of the data provided
(normally a subset).

`complexity()`: Returns a `vector` with the binary split indexes.

## Examples

``` r
normalized <- znorm(motifs_discords_small)
fake_data <- c(rep(3, 100), rep(2, 100), rep(1, 100))
correlation <- ed_corr(fake_data, 50)
fake_data <- c(rep(0.5, 100), rep(1, 100), rep(0.1, 100))
euclidean <- corr_ed(fake_data, 50)
fake_data <- c(1, 1, 4, 5, 2, 3, 1, 7, 9, 4, 5, 2, 3)
mode <- mode(fake_data)
fake_data <- c(1, 1.4, 4.3, 5.1, 2, 3.6, 1.24, 2, 9, 4.3, 5, 2.1, 3)
res <- std(fake_data)
fake_data <- c(1, 1.4, 4.3, 5.1, 2, 3.6, 1.24, 1, 9, 4.3, 5, 2.1, 3)
res <- normalize(fake_data)
fake_data <- c(1, 1.4, 4.3, 5.1, 2, 3.6, 1.24, 8, 9, 4.3, 5, 2.1, 3)
res <- complexity(fake_data)
fake_data <- c(10)
res <- binary_split(fake_data)
```
