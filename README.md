TFBoost: a package for a tree-based functional boosting algorithm
================
Xiaomeng Ju and Matias Salibian Barrera
2021-06-15

This repository contains `R` code implementing a tree-based boosting
algorithm for scalar-on-function regression. The code provides a fit for
a multi-index regression model introduced in (add reference)…

## Install and load package

You can install the development version of the package in R using:

``` r
devtools::install_github("xmengju/TFBoost", auth_token = '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')
```

Once installed you can load the package with:

``` r
library(TFBoost)
```

## Example: German electricity data

Below we illustrate the use of the package with the German electricity
dataset. The original data is provided in the on-line supplementary
materials of [Liebl
(2013)](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-7/issue-3/Modeling-and-forecasting-electricity-spot-prices--A-functional-data/10.1214/13-AOAS652.full).The
data consist of electricity spot prices traded at the European Energy
Exchange (EEX) in Leipzig and electricity demand reported by European
Network of Transmission System Operators for Electricity from January
1st 2006 to September 30th 2008. We removed weekends and holidays from
the data and ended up with data collected on 638 days and provided it in
the package as `GED`. We treated hourly evaluated electricity spot
prices (`GED$price`) as the predictor, represented as a vector of
dimension 24 for each day, and the daily average of electricity demand
(`GED$demand`) as the response.

We load the data and plot the first 10 predictor curves.

``` r
data(GED)
matplot(t(GED$price[1:10,]), lty = 1, type = "l", ylab = "price", xlab = "hour")
```

![](README_files/figure-gfm/plot-1.png)<!-- -->

In order to train our predictor, we split the data set into a `training`
set (with 60% of the available data), a `validation` set and a `test`
set (both with 20% of the data). We first randomly select the
observations for each of these three sets:

``` r
n <- nrow(GED$price) 
n0 <- floor( 0.2 * n) 
set.seed(123)
idx_test <- sample(n, n0)
idx_train <- sample((1:n)[-idx_test], floor( 0.6 * n ) )
idx_val <- (1:n)[ -c(idx_test, idx_train) ] 
```

We now create the matrices of explanatory variables (`x`) and vectors of
responses (`y`) corresponding to this partition.
<!-- Note that `ytrain` and `yval` may contain outliers. -->

``` r
xtrain <- GED$price[idx_train, ]
ytrain <- GED$demand[idx_train ]
xval <- GED$price[idx_val, ]
yval <- GED$demand[idx_val ]
xtest <- GED$price[idx_test, ]
ytest <- GED$demand[idx_test ]
```

The `TFBoost` function implements the proposed tree-based functional
boosting algorithm with two options for the base learner: type A tree or
type B tree.

We now explain how to fit a `TFBoost` estimator and compare it with the
`fgam` estimator proposed in [McLean el
al. (2014)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3982924/) and
implemented in the `refund` package. To specify the tree type, the user
needs to set `tree_type = A` or `tree_type = B` in `tree_control`.
Below, we will fit `TFBoost` with the type B tree, which trains much
faster compared to the type A tree.



The depth of the base learners in `TFBoost` is set with the argument
`max_depth`. The argument `niter` specifies the number of iterations to
be used (we set `niter = 1000`).
