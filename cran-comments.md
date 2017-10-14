## Test environments
* local OS X install, R 3.4.2
* ubuntu 14.04 (on travis-ci)
  * R 3.3.3 (oldrel)
  * R 3.4.2 (release)
  * R Under development (unstable) (2017-10-13 r73550) (devel)
* win-builder (devel and release)

## R CMD check results
* There were no ERRORs or WARNINGs.
* There were no NOTES on the local (Mac) or travis-ci builds.
* However, there was 1 NOTE from win-builder -- both words are surnames and are correct:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'John A. Ramey <johnramey@gmail.com>'

Possibly mis-spelled words in DESCRIPTION:
  Dudoit (11:16)
  Ramey (10:10)
  al (10:19, 11:26, 12:5, 12:68)
  et (10:16, 11:23, 11:76, 12:65)

## Reverse dependencies

* The `mlr` package suggests `sparsediscrim`. There were no errors or warnings.
