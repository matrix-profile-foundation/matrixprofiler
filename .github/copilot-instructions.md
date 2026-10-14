# AI Coding Agent Instructions for matrixprofiler

**matrixprofiler** is an R package implementing the Matrix Profile algorithms for time series analysis, a core toolkit for motif discovery, anomaly detection, and time series similarity. It separates performance-critical computation (C++ in `src/`) from R wrapper logic.

## Persistent Memory

Whenever the user indicates they are starting work for the day, consult the `.github/copilot-memory.md` file to obtain the most recent information about the project status.
When the user indicates they have finished work for the day, update the `.github/copilot-memory.md` file with the most recent information about the project status.

## Architecture & Key Components

### R/C++ Integration Pattern
- **R wrappers** ([R/*.R](R/)) expose user-facing functions with extensive parameter validation
- **C++ implementations** ([src/*.cpp](src/)) contain compute-intensive algorithms via Rcpp
- Critical: Each R function has corresponding `*_rcpp()` export (e.g., `mpx()` calls `mpx_rcpp()` or `mpx_rcpp_parallel()`)
- All new algorithms need both R and C++ versions

### Parallel Processing Architecture
- Uses **RcppParallel** with TBB backend (default, can switch to tinythread)
- Backend selection happens in [R/zzz.R](.onAttach/zzz.R) via `RCPP_PARALLEL_BACKEND`
- Single vs. multi-threaded dispatch in R layer: `if (n_workers > 1)` calls `*_parallel()` variant
- Thread management: `RcppParallel::setThreadOptions(numThreads = n_workers)` wraps parallel calls (see [R/mpx.R](R/mpx.R#L82-L91))
- C++ parallel workers: `struct Worker : public RcppParallel::Worker` pattern with `operator()(begin, end)` override
- Conditional compilation: `#if RCPP_PARALLEL_USE_TBB` for mutex types (tbb::spin_mutex vs tthread::mutex) in [src/mpx_parallel.cpp](src/mpx_parallel.cpp#L24)

### Core Algorithms (STAMP, STOMP, SCRIMP, MPX)
Each has:
- **Single-threaded**: `*_rcpp()` - base computation
- **Multi-threaded**: `*_rcpp_parallel()` - RcppParallel Worker pattern
- **AB-join variants**: For cross-sequence similarity (e.g., `mpxab_rcpp_parallel()`)
- **Pre-computation functions**: MASS variants use precomputed statistics to optimize distance profiles

## Code Conventions & Patterns

### R Function Structure
All user functions follow this template (see [R/mass.R](R/mass.R#L27-L77)):
1. Parse & validate arguments with `checkmate::qassert()` - includes debug statements `"!!!DEBUG Parsing Arguments"`
2. Initialize `result <- NULL` and register `on.exit()` cleanup that returns result if not null
3. Wrap computation in `tryCatch()` to catch errors gracefully
4. Return via `on.exit()` exit handler so cleanup happens automatically

**Example error-handling pattern:**
```r
result <- NULL
on.exit({
  if (is.null(result)) return(invisible(NULL))
  else { result$ez <- ez; return(result) }
}, TRUE)
tryCatch({ result <- function_call() }, error = print)
```

### C++ Patterns
- **Window functions**: Precompute stats (mean, std, fft) before distance calculations
- **Progress bars**: RcppProgress Progress(100, progress) across parallelFor loops
- **User interruption**: Catch `RcppThread::UserInterruptException` to allow graceful stopping with partial results
- **Distance conversion**: Normalize/denormalize between correlation and Euclidean metrics (see [src/mpx_parallel.cpp](src/mpx_parallel.cpp#L447-L451))

### Data Input Validation
- Coerce to numeric: `data <- as.numeric(data)` then `checkmate::qassert(data, "N+")` (N+ = numeric, positive)
- Window size validation: Time series must be > 2× window_size
- Exclusion zone: Prevents trivial matches; default 0.5 × window_size

## Testing Strategy

### Test Organization
- [tests/testthat/](tests/testthat/) contains **snapshot tests** with reference data
- `motifs_discords_small` is test dataset loaded from [data/motifs_discords_small.rda](data/motifs_discords_small.rda)
- Tests use `expect_snapshot_value()` for numerical outputs to catch regressions
- Parallel vs. serial test parity: Tests run algorithms with `n_workers=1` and `n_workers=2`, assert results match with `expect_equal()`

### Running Tests
- Use VS Code task "R: test units": `testthat::test_local(reporter = testthat::SummaryReporter)`
- Linting: "R: lint package" (styler) + "R: lint all code" (clang-format for C++)
- Code analysis: "PVS-Analyze" and "PVS-Convert" tasks for static analysis

## Development Workflow

### Build & Development
- **Build package**: "R: build current package" task runs `devtools::load_all()`
- **C++ compilation**: Happens automatically via Rcpp; to rebuild just one file use clang-format
- **Makevars**: [src/Makevars](src/Makevars) enables TBB + RcppParallel linkage

### Debugging & Validation
- Enable debugme: [R/zzz.R](R/zzz.R#L2) calls `debugme::debugme()` on load
- Use "!!!DEBUG" strings (exclamation marks enable them in debugme) for development logging
- ASAN/UBSAN checks: "R: ASAN/USBAN" task detects memory errors in C++ (Windows-oriented but useful)

### Documentation
- Function docs: roxygen2 comments in R files (`#' @param`, `#' @details`, `#' @examples`)
- Package docs generated via roxygen2; rebuild with VS Code Rcpp extension
- Versioning: Update [DESCRIPTION](DESCRIPTION) before release

## Common Tasks

### Adding a New Algorithm
1. Create R wrapper in [R/newname.R](R/) following mass.R/mpx.R template
2. Implement C++ in [src/newname.cpp](src/) and [src/newname.h](src/)
3. Create both serial `newname_rcpp()` and `newname_rcpp_parallel()` exports
4. Add RcppExports via roxygen: `#' @useDynLib matrixprofiler, .registration = TRUE`
5. Add test in [tests/testthat/test-newname.R](tests/testthat/) comparing single vs. parallel
6. Use VS Code task "R: build current package" to recompile

### Fixing Parallel Issues
- Check thread backend in [R/zzz.R](R/zzz.R#L10): default TBB, set `RCPP_PARALLEL_BACKEND` env var to switch
- Verify `setThreadOptions()` wraps parallel call in R code
- Ensure C++ Worker struct inherits from `RcppParallel::Worker` and implements `operator()(begin, end)`
- Use `RcppParallel::parallelFor()` for TBB or `RcppParallel2::ttParallelFor()` for tinythread (conditional via macro)

### Performance Optimization
- Profile with proffer package (listed in dev/make.R)
- Pre-compute expensive operations (statistics, FFTs) before parallelFor loops
- Grain size matters: 4×window_size typical (see [src/mpx_parallel.cpp](src/mpx_parallel.cpp#L165))
- Use `RVector<T>` accessor in Worker structs for efficient thread-safe reads

## Key Files Reference
- [R/matrixprofiler-package.R](R/matrixprofiler-package.R) - package overview & parallel backend docs
- [src/mpx_parallel.cpp](src/mpx_parallel.cpp) - exemplary Worker pattern with TBB/tinythread handling
- [R/mpx.R](R/mpx.R) - exemplary R wrapper with parallel dispatch
- [tests/testthat/test-stamp.R](tests/testthat/test-stamp.R) - exemplary test showing parity checks
