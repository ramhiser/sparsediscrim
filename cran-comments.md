## Test environments
* local OS X install, R 3.4.1
* ubuntu 12.04 (on travis-ci)
  * R 3.3.3 (oldrel)
  * R 3.4.1 (release)
  * R Under development (unstable) (2017-08-11 r73085) (devel)
* win-builder (devel and release)

## R CMD check results
* There were no ERRORs or WARNINGs.
* There were no NOTES on the local (Mac) or travis-ci builds.
* However, there was 1 NOTE from win-builder:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'John A. Ramey <johnramey@gmail.com>'

New submission

Package was archived on CRAN

CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2017-04-21 as check errors were not
    corrected despite reminders.

## Reverse dependencies

* The `mlr` package suggests `sparsediscrim`. There were no errors or warnings.
