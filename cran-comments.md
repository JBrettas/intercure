## Resubmission
This is a resubmission. In this version I have:

* Used parallel and doParallel instead of snow and doSNOW
* Controlled the NAMESPACE using stats4::vcov and stats4::coef
* Reduced the examples runtime
* Removed some unnecessary unit tests

## Test environments
* local windows 8.1 pro, R 3.2.2
* ubuntu 12.04 (on travis-ci), R 3.2.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Julio Brettas <jbrettas@ime.usp.br>'
New submission

## Downstream dependencies
The package is new and have no downstream dependencies.
