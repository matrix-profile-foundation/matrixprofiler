# CRAN Comments

## Comments for this version (0.1.10)

C++11 flags are deprecated and defaults to C++17, so this requirement was removed from the package.
Also, some old links were fixed.

## Test environments

- rhub: all platforms (\*nix, solaris, macos, windows, with gcc and clang), all ok,
  -  3 [VM] macos          R-* (any version)                     macos-13 on GitHub
  -  9 [CT] clang-ubsan    R-devel (2025-11-30 r89082)           Ubuntu 22.04.5 LTS
  - 10 [CT] clang16        R-devel (2025-11-29 r89077)           Ubuntu 22.04.5 LTS
  - 16 [CT] gcc-asan       R-devel (2025-11-30 r89082)           Fedora Linux 40 (Container Image)
  - 27 [CT] ubuntu-clang   R-devel (2025-11-30 r89082)           Ubuntu 22.04.5 LTS
  - 30 [CT] ubuntu-release R-4.5.2 (2025-10-31)                  Ubuntu 24.04.3 LTS
  - 31 [CT] valgrind       R-devel (2025-11-30 r89082)           Fedora Linux 38 (Container Image)

- win-builder: devel, release, oldrel

## R CMD check results

── R CMD check results ──────────────────────────────────────── matrixprofiler 0.1.10 ────
Duration: 1m 17s

❯ checking installed package size ... NOTE
    installed size is 13.1Mb
    sub-directories of 1Mb or more:
      libs  12.9Mb

❯ checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

❯ checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-Wdate-time’ ‘-Werror=format-security’ ‘-Wformat’

0 errors ✔ | 0 warnings ✔ | 3 notes ✖

## Downstream dependencies

- No reverse dependencies yet

## Known Issues (a.k.a NOTES)

- GNU make is a SystemRequirements.

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
