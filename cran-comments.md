## Comments

New package implementing the core functions of `tsmp` package.
That package will then imports this as a dependency.

## Test environments

- GitHub Actions (ubuntu-16.04): devel, release, oldrel, 3.5, 3.4, 3.3
- GitHub Actions (windows): release, oldrel
- GitHub Actions (macOS): release
- win-builder: devel

## R CMD check results

`0 errors | 0 warnings | 0 notes`

## Downstream dependencies

- No reverse dependencies yet

## Known Issues (a.k.a NOTES)

- Found the following (possibly) invalid file URI:
  URI: .github/CODE_OF_CONDUCT.md
  From: README.md

  - This is ok.

- GNU make is a SystemRequirements.

  - Requirement of package RcppParallel. I haven't find a workaround to solve this NOTE.

- Installed size is 7.5Mb.

  - This is due to datasets in this package. I believe they are essential to learning all the features
    of this package.

- Uses the superseded package: `doSNOW`

  - `doSNOW` has a property that allows to use progress bar that `parallel` does not.
  - Working in finding a better solution to drop this dependency. Not found yet.

- (possibly) invalid URLs: https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html

  - Debian: libcurl throws an error on certificate check. Nothing to do about this.

- Authors@R field gives persons with non-standard roles
  - These non-standard roles where appropriately chosen using [MARC Code List for Relators](https://www.loc.gov/marc/relators/relaterm.html)
