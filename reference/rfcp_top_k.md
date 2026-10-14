# Top-K Relative Frequency Contrast Profile

Discovers the top-K contrasting patterns (Platos) between a positive
time series and a negative time series using an iterative Relative
Frequency Contrast Profile (RFCP) procedure.

## Usage

``` r
rfcp_top_k(
  positive_data,
  negative_data,
  window_size,
  k_platos = 2L,
  max_freq = 5L,
  exclusion_radius = NULL,
  min_rms = 1e-04,
  adaptive_max_freq = FALSE,
  min_freq = 2L,
  n_workers = 1L,
  backend = "cpu",
  progress = FALSE,
  query_block_size = NULL,
  progress_callback = NULL,
  checkpoint_callback = NULL,
  resume_state = NULL
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

- k_platos:

  Integer number of distinct Platos to discover. Defaults to `2L`.

- max_freq:

  Integer scalar or integer vector specifying the number of
  non-overlapping neighbor ranks requested for each Plato iteration. If
  a scalar is provided, it is used for all iterations. If a vector is
  provided, its length must equal `k_platos`, allowing dynamic
  re-evaluation of frequency for rarer morphologies. Defaults to `5L`.

- exclusion_radius:

  Integer exclusion radius between accepted neighbor starts and for
  padding concatenated segments. `NULL` uses `window_size`.

- min_rms:

  Numeric threshold for root-mean-square contrast. Iteration stops early
  if the detected Plato has an RMS contrast below this value. Defaults
  to `1e-4`.

- adaptive_max_freq:

  Logical. If `TRUE`, dynamically adapts `max_freq` for subsequent
  iterations based on the observed support (number of valid positive
  self-join neighbors with distance below \\\sqrt{2 \cdot
  \text{window\\size}}\\) of the current Plato. Defaults to `FALSE`.

- min_freq:

  Integer minimum rank frequency floor when `adaptive_max_freq = TRUE`.
  Defaults to `2L`.

- n_workers:

  Number of native worker threads for the CPU backend. Ignored when
  `backend = "cuda"`. Defaults to `1L`. For parallel execution across
  processes, prefer a PSOCK cluster and use `n_workers = 1` within each
  worker. See
  [matrixprofiler](https://matrix-profile-foundation.github.io/matrixprofiler/reference/matrixprofiler.md)
  for fork limitations.

- backend:

  Computation backend: `"cpu"` (default) or `"cuda"`.

- progress:

  Logical progress flag for each RFCP iteration. Defaults to `FALSE`.

- query_block_size:

  Optional integer query block size for chunked execution. If `NULL`
  (default), queries are processed monolithically in a single call. If
  specified, queries are processed in contiguous blocks of size
  `query_block_size`.

- progress_callback:

  Optional callback function `function(event)` invoked periodically with
  structured progress events (`"round_start"`, `"block_done"`,
  `"round_done"`).

- checkpoint_callback:

  Optional callback function `function(state)` invoked after each block
  and Plato with a serializable state list.

- resume_state:

  Optional state list from a previous execution (`checkpoint_callback`)
  used to resume execution idempotently.

## Value

An object of class `"rfcp_top_k"`, which is a list containing:

- `k_requested`: requested number of Platos.

- `k_found`: number of Platos discovered before stopping.

- `window_size`: subsequence length.

- `exclusion_radius`: exclusion radius used.

- `adaptive_max_freq`: whether adaptive frequency was enabled.

- `min_freq`: minimum frequency floor used.

- `max_freq`: vector of frequencies used for each Plato.

- `platos`: list of discovered Plato objects, each containing:

  - `k`: Plato iteration rank (1, 2, ...).

  - `plato_index`: one-based start index in `positive_data`.

  - `plato_rms`: contrast score at `plato_index`.

  - `plato_subsequence`: numeric vector containing the pattern.

  - `plato_twin_index`: one-based start index of nearest self-join
    neighbor.

  - `rfcp_at_plato`: clipped contrasts across accepted ranks.

  - `rfmp_aa_at_plato`: positive self-join distances.

  - `rfmp_ab_at_plato`: negative AB distances.

  - `rfmpi_aa_at_plato`: positive neighbor indices.

  - `rfmpi_ab_at_plato`: negative neighbor indices.

  - `effective_max_freq`: highest paired rank observed.

- `negative_data_final`: the negative series augmented with the
  suppressed Plato occurrences and exclusion zones.

## Details

In each iteration \\i \in \\1, \dots, K\\\\, the algorithm:

1.  Executes
    [`rfcp`](https://matrix-profile-foundation.github.io/matrixprofiler/reference/rfcp.md)
    using the current augmented negative series and the requested
    frequency `max_freq[i]`.

2.  Identifies the highest-contrast query window (\\\text{Plato}\_i\\)
    outside the exclusion zone of previously discovered Platos and their
    associated neighbors.

3.  If the contrast score is below `min_rms` or no valid candidate
    remains, iteration terminates early.

4.  Extracts the subsequences of \\\text{Plato}\_i\\ and all its
    accepted positive self-join neighbors (`rfmpi_aa_at_plato`),
    expanded by the exclusion zone: \$\$\[\max(1, j -
    \text{exclusion\\radius}), \min(N^{(+)}, j + \text{window\\size} -
    1 + \text{exclusion\\radius})\]\$\$

5.  Merges overlapping intervals to maintain signal continuity, and
    appends the resulting segments to the negative series separated by
    `NA` barriers. Because `rfcp` is NA-aware, this suppresses the
    contrast of that morphology and its trivial overlaps in subsequent
    iterations without introducing spurious cross-boundary windows.

6.  If `adaptive_max_freq = TRUE`, re-evaluates the frequency for the
    next iteration based on the count of finite AA neighbors with
    distance below \\\sqrt{2 \cdot \text{window\\size}}\\, bounded below
    by `min_freq`.

## Examples

``` r
set.seed(42)
# Create positive series with 5 sine pulses and 2 cosine pulses
t <- seq(0, 4 * pi, length.out = 20)
sine_pulse <- sin(t)
cos_pulse <- cos(t)
pos <- rnorm(300, sd = 0.05)
for (p in c(20, 60, 100, 140, 180)) pos[p:(p + 19)] <- pos[p:(p + 19)] + sine_pulse
for (p in c(220, 260)) pos[p:(p + 19)] <- pos[p:(p + 19)] + cos_pulse * 1.5
neg <- rnorm(300, sd = 0.05)

top_platos <- rfcp_top_k(pos, neg, window_size = 20, k_platos = 2, max_freq = c(4, 2))
print(top_platos)
#> Top-2 Relative Frequency Contrast Profile (RFCP)
#> Window size: 20 | Exclusion radius: 20 | Requested: 2 | Adaptive: FALSE
#>   Plato 1: index = 97, RMS = 0.683502, max_freq = 4 (twin = 177)
#>   Plato 2: index = 279, RMS = 0.341239, max_freq = 2 (twin = 239)
```
