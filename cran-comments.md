## Test environments
* local OS X install, R 3.4.0
* ubuntu 12.04 (on travis-ci)
  * R 3.2.5 (oldrel)
  * R 3.3.0 (release)
  * R 3.4.0 (devel)
* win-builder (devel and release)

## R CMD check results
* There were no ERRORs or WARNINGs.
* There were no NOTES on the local (Mac) or travis-ci builds.
* However, there was 1 NOTE from win-builder:

Maintainer: 'John A. Ramey <johnramey@gmail.com>'

License components with restrictions and base license permitting such:
MIT + file LICENSE
File 'LICENSE':
YEAR: 2016
  COPYRIGHT HOLDER: John A. Ramey <johnramey@gmail.com>

## Reverse dependencies

* The `mlr` package suggests `sparsediscrim`. There were no errors or warnings.
