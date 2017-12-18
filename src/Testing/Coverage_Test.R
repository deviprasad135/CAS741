library(devtools)
library(covr)
library(DT)
library(testthat)
devtools::test()                     #unit testing
cov <- package_coverage()            #code coverage
report(cov)                          #summary for code coverage
