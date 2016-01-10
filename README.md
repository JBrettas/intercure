<!-- README.md is generated from README.Rmd. Please edit that file -->
intercure
=========

[![Build Status](https://travis-ci.org/JBrettas/intercure.svg?branch=master)](https://travis-ci.org/JBrettas/intercure)

The intercure package provides R implementations for different cure rate models to interval censored data. The goal is to estimate the effects of each given covariate using different specifications and computing the estimated cure fraction with it. For simulation and illustrative examples, the package also provides functions to generate a dataset using one of the cure rate specification.

Two models for interval censored data are considered on this package: the bounded cumulative hazard model proposed on Shen and Hao (2009); the frailty model proposed on Lam, Wong, and Zhou (2013), with its extension to clustered data presented on Lam and Wong (2014).

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

References
----------

Lam, Kwok Fai, and Kin Yau Wong. 2014. “Semiparametric Analysis of Clustered Interval-Censored Survival Data with a Cure Fraction.” *Computational Statistics and Data Analysis* 79: 165–74.

Lam, Kwok Fai, Kin Yau Wong, and Feifei Zhou. 2013. “A Semiparametric Cure Model for Interval-Censored Data.” *Biometrical Journal* 55 (5): 771–88.

Shen, Yu, and Liu Hao. 2009. “A Semiparametric Regression Cure Model for Interval-Censored Data.” *Journal of the American Statistical Association* 104 (487): 1168–78.
