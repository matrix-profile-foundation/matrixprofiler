## Comments

New package implementing the core functions of 'tsmp' package.
That package will then import this as a dependency.

## Test environments

- GitHub Actions (ubuntu-16.04): devel, release, oldrel
- GitHub Actions (windows): release, oldrel
- GitHub Actions (macos): release
- rhub: all platforms (*nix, solaris, macos, windows, with gcc and clang), all ok,
  except 'solaris-x86-patched-ods' in which a dependency fails, not this package.
- Raspbian (armv7l-linux-gnueabihf, 32-bits): release
- win-builder: devel, release, oldrel

## R CMD check results

`0 errors | 0 warnings | 1 notes`

## Downstream dependencies

- No reverse dependencies yet

## Known Issues (a.k.a NOTES)

- GNU make is a SystemRequirements.

## Last CRAN comments

- Last known issue:
  "runtime error: -inf is outside the range of representable values of type 'int'"
  
  - Hope it is fixed.


## Old comments

- Please use only undirected quotation marks in the description text.
  e.g. `tsmp` --> 'tsmp' (fixed)

- Please add \value to .Rd files regarding exported methods and explain
  the functions results in the documentation. Please write about the
  structure of the output (class) and also what the output means. (If a
  function does not return a value, please document that too, e.g.
  \value{No return value, called for side effects} or similar) (fixed, all functions has returns value explained)

- Missing Rd-tags:
         pipe.Rd: \value (dropped, will be keep in the 'tsmp' package, not this one)
         test_equal.Rd: \value (dropped, internal use only)

- Please always add all authors, contributors and copyright holders in the
  Authors@R field with the appropriate roles.
  e.g.: Jarek Tuszynski (read below)

  - Authors@R field gives persons with non-standard roles
  - These non-standard roles where appropriately chosen using [MARC Code List for Relators](https://www.loc.gov/marc/relators/relaterm.html)

