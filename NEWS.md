NEWS
================
Francisco Bischoff

25 Jan 2023

<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->

## matrixprofiler 0.1.9

-   Small fixes for cran

## matrixprofiler 0.1.8

-   Fixed memory leak on mpx_rcpp_parallel().
-   Added contrast profile function: contrast().
-   Several code lint.

## matrixprofiler 0.1.7

-   Fixed issue on Makefile for Unix/Mac systems that triggered some
    CRAN policy about “Packages should not attempt to disable compiler
    diagnostics, nor to remove other diagnostic information such as
    symbols in shared objects.” The objective was to reduce the compiled
    library that even in release mode have being compiled in debug mode.
    More info at <http://dirk.eddelbuettel.com/blog/2017/08/14/>

## matrixprofiler 0.1.5

-   Fixed `Rcpp` as specified by the PR \#15 from Dirk @eddelbuettel

## matrixprofiler 0.1.4

-   Dropped PAA algorithm
-   Added some math helper functions

## matrixprofiler 0.1.3

-   CRAN fixes

## matrixprofiler 0.1.0

-   Initial implementation of this package.
-   This package will keep all core functions that will allow you to use
    the Matrix Profile concept as a toolkit.
-   It will be the main dependency of the already available package
    `tsmp`.
