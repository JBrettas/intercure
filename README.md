<!-- README.md is generated from README.Rmd. Please edit that file -->
intercure
=========

[![Build Status](https://travis-ci.org/JBrettas/intercure.svg?branch=master)](https://travis-ci.org/JBrettas/intercure)

The intercure package provides R implementations for different cure rate models to interval censored data. The goal is to estimate the effects of each given covariate using different specifications and computing the estimated cure fraction with it. For simulation and illustrative examples, the package also provides functions to generate a dataset using one of the cure rate specification.

Two models for interval censored data are considered on this package: the bounded cumulative hazard model proposed on Liu and Shen (2009); the frailty model proposed on Lam, Wong, and Zhou (2013), with its extension to clustered data presented on Lam and Wong (2014).

Installing
----------

``` r
devtools::install_github("JBrettas/intercure")
```

Functions
---------

To fit a cure rate model:

-   inter\_bch
-   inter\_frailty
-   inter\_frailty\_cl

To generate datasets:

-   sim\_bch
-   sim\_frailty
-   sim\_frailty\_cl

How to use
----------

Click [here](http://rpubs.com/JBrettas/howtointercure) for a tutorial on how to use the package.

Special thanks
--------------

My greatest thanks goes to Professor Lam, who gave me help on theorical details of the frailty model and its implementation. Also, my special thanks to Professor Shen, for all the answered doubts and for sharing the original C routine. And of course to Professor Tunes, who greatly oriented me on my master's degree.

References
----------

Lam, Kwok Fai, and Kin Yau Wong. 2014. “Semiparametric Analysis of Clustered Interval-Censored Survival Data with a Cure Fraction.” *Computational Statistics and Data Analysis* 79: 165–74.

Lam, Kwok Fai, Kin Yau Wong, and Feifei Zhou. 2013. “A Semiparametric Cure Model for Interval-Censored Data.” *Biometrical Journal* 55 (5): 771–88.

Liu, Hao, and Yu Shen. 2009. “A Semiparametric Regression Cure Model for Interval-Censored Data.” *Journal of the American Statistical Association* 104 (487): 1168–78.
