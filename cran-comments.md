# CRAN Comments

## Comments for this version (0.1.8)

Due to CRAN policy, I had to remove the Makefile command --strip-debug.
Therefore, the library now has more than 15MB in size on *nix systems.

Removed VignetteBuilder field since this package does not contain a vignette.

## Test environments

-   GitHub Actions (ubuntu-16.04): devel, release, oldrel
-   GitHub Actions (windows): release, oldrel
-   GitHub Actions (macos): release
-   rhub: all platforms (*nix, solaris, macos, windows, with gcc and clang), all ok,
  except:
    -   `solaris-x86-patched-ods` in which a dependency fails, not this package.
    -   `linux-x86_64-centos-epel`: not available: 'spelling'; found 'abort'; found 'printf'. I could not find why the code behaves like this on this distribution.
    -   `ubuntu-rchk`: Bioconductor does not yet build and check packages for R version 4.2; see <https://bioconductor.org/install>. There is nothing in the package code that asks for R 4.2 version.
-   Raspbian (armv7l-linux-gnueabihf, 32-bits): release
-   win-builder: devel, release, oldrel

## R CMD check results

── R CMD check results ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── matrixprofiler 0.1.7.9000 ────
Duration: 2m 15.1s

❯ checking installed package size ... NOTE
    installed size is 16.5Mb
    sub-directories of 1Mb or more:
      libs  16.3Mb

❯ checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

## Downstream dependencies

-   No reverse dependencies yet

## Known Issues (a.k.a NOTES)

-   GNU make is a SystemRequirements.

## Current comments

-   None yet.

## Old comments

-   [x] I see that last version was archived due to OSX bug on compiling cpp code. I wasn't warned about it, at least I didn't see any e-mail in my inbox. Nevertheless, the "bug" is known and I fixed it in this version. Already tested on rhubs.

-   [x] Please use only undirected quotation marks in the description text.
  e.g. `tsmp` --> 'tsmp'

-   [x] Please add \value to .Rd files regarding exported methods and explain
  the functions results in the documentation. Please write about the
  structure of the output (class) and also what the output means. (If a
  function does not return a value, please document that too, e.g.
  \value{No return value, called for side effects} or similar) (fixed, all functions has returns value explained)

-   [x] Missing Rd-tags:
         pipe.Rd: \value (dropped, will be keep in the 'tsmp' package, not this one)
         test_equal.Rd: \value (dropped, internal use only)

-   [x] Please always add all authors, contributors and copyright holders in the
  Authors@R field with the appropriate roles.
  e.g.: Jarek Tuszynski (read below)

    -   Authors@R field gives persons with non-standard roles
    -   These non-standard roles where appropriately chosen using [MARC Code List for Relators](https://www.loc.gov/marc/relators/relaterm.html)
