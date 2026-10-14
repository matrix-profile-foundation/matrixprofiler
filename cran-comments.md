# CRAN Comments

## Comments for this version (0.1.12)

- Complete migration from Intel TBB / `RcppParallel` to `RcppThread`: all multithreaded algorithms (`scrimp`, `stomp`, `mpx`, `mstomp`, `mass`) are now cleanly threaded with `RcppThread` without external TBB linkage issues or fork-unsafe crashes.
- Implemented high-performance native C++ `mSTOMP` engine (`mstomp()`), yielding dramatic speedups for multidimensional matrix profile computation.
- Rolling window functions (`mov_mean`, `mov_sum`, `mov_var`, `mov_std`, `mov_min`, `mov_max`, `movmean_std`, `muinvn`) now support explicit `na_strategy = c("window_na", "remove", "error")` and standard `na.rm`.
- Extended multithreading to all MASS variants (`v2`, `weighted`, `non_normalized`, `absolute`).
- Standardized invalid profile indices to `NA_integer_` across algorithms.
- `R CMD check --as-cran` passed with 0 errors and 0 warnings.
- "The NOTE ‘-mno-omit-leaf-frame-pointer’ appears only on local Ubuntu build environment, not on CRAN test machines."

## Test environments

- rhub: linux-devel, windows-devel, macos-devel, clang-ubsan, valgrind

- win-builder: devel, release

## R CMD check results

── R CMD check results ───────────────────────────────────────── matrixprofiler 0.1.12 ────
Duration: 1m 18.9s

❯ checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-mno-omit-leaf-frame-pointer’

0 errors ✔ | 0 warnings ✔ | 1 note ✖

## Downstream dependencies

- No reverse dependencies yet

## Known Issues (a.k.a NOTES)

- GNU make is a SystemRequirements.
- Package size

## Current comments

- None yet.

## Old comments

## Comments for this version (0.1.9)

- [x] Due to CRAN policy, I had to remove the Makefile command --strip-debug. Therefore, the library now has more than
      15MB in size on \*nix systems. Removed VignetteBuilder field since this package does not contain a vignette.

- [x] I see that last version was archived due to OSX bug on compiling cpp code. I wasn't warned about it, at least I
      didn't see any e-mail in my inbox. Nevertheless, the "bug" is known and I fixed it in this version. Already tested on
      rhubs.

- [x] Please use only undirected quotation marks in the description text.
      e.g. `tsmp` --> 'tsmp'

- [x] Please add \value to .Rd files regarding exported methods and explain
      the functions results in the documentation. Please write about the
      structure of the output (class) and also what the output means. (If a
      function does not return a value, please document that too, e.g.
      \value{No return value, called for side effects} or similar) (fixed, all functions has returns value explained)

- [x] Missing Rd-tags:
      pipe.Rd: \value (dropped, will be keep in the 'tsmp' package, not this one)
      test_equal.Rd: \value (dropped, internal use only)

- [x] Please always add all authors, contributors and copyright holders in the
      Authors@R field with the appropriate roles.
      e.g.: Jarek Tuszynski (read below)

  - Authors@R field gives persons with non-standard roles
  - These non-standard roles where appropriately chosen using
    [MARC Code List for Relators](https://www.loc.gov/marc/relators/relaterm.html)
