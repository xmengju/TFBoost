TFBoost: a package for a tree-based functional boosting algorithm
================
Xiaomeng Ju and Matias Salibian Barrera
2021-06-17

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

## An example: German electricity data

Below we illustrate the use of the package with the German electricity
dataset. The original data is provided in the on-line supplementary
materials of [Liebl
(2013)](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-7/issue-3/Modeling-and-forecasting-electricity-spot-prices--A-functional-data/10.1214/13-AOAS652.full).The
data consist of electricity spot prices traded at the European Energy
Exchange (EEX) in Leipzig and electricity demand reported by European
Network of Transmission System Operators for Electricity from January
1st 2006 to September 30th 2008. We excluded weekends and holidays from
the data and ended up with data collected on 638 days and provided it in
the package as `GED`. We treated hourly evaluated electricity spot
prices (`GED$price`) as the predictor, represented as a vector of
dimension 24 for each day, and the daily average of electricity demand
(`GED$demand`) as the response.

We plot the data, highlighting 10 randomly chosen predictor curves in
red. We also plot the average daily demand.

``` r
data(GED)
matplot(t(GED$price), lty = 1, type = "l", ylab = "price", xlab = "hour", 
        main = "hourly electricity price (10 days)", col='gray70', ylim=c(0, 1000))
set.seed(123)
n <- nrow(GED$price)
matplot(t(GED$price[sample(n, 10), ]), type='l', col='red', add=TRUE, lty=1)
plot(GED$date,GED$demand, type = "l", ylab = "demand", xlab = "date", 
     main = "average daily electricity demand")
```

<img src="README_files/figure-gfm/plot-1.png" width="45%" /><img src="README_files/figure-gfm/plot-2.png" width="45%" />

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
al. (2014)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3982924/) and
implemented in the `refund` package. To specify the tree type, the user
needs to set `tree_type = A` or `tree_type = B` in `tree_control`.
Below, we will fit `TFBoost` with the type B tree, which trains much
faster compared to the type A tree.

The following are parameters required for our estimator

``` r
tree_type  <- "B" # type of the base learner
num_dir <- 20  # number of random directions for type B tree
gg <- 1:24  # specify the grid the functional predictor was evaluated on
tt <- c(0,24) # domain of the functional predictor
niter <- 1000 # number of boosting iterations 
make_prediction <- TRUE # make predictions based on test data
loss <-  "l2" # loss for the boosting algorithm ("l2", "lad", or one specified by user_func)
shrinkage <- 0.05 # shrinkage parameter for boosting
nknot <- 3 # the number of interior knots for cubic B-spline basis
```

The depth of the base learners in `TFBoost` is set with the argument
`d`. We considered `d` from 1 to 4, and chose the depth the minimizes
the mean-squared-error on the `validation set` at early stopping time.
Below we train `TFBoost` and select the depth `d`. This step may take
several minutes to run:

``` r
tree.depths <- 1:4
model_TFBoost_list <- vector('list', length(tree.depths))
val_error_vec <- rep(NA, length(tree.depths))
for(dd in 1:length(tree.depths)) {
  model_TFBoost_list[[dd]] <- 
    TFBoost(x_train = xtrain, y_train = ytrain,  x_val = xval,  y_val = yval,
            x_test = xtest, y_test = ytest, grid = gg, t_range  = tt, niter = niter,
            control = TFBoost.control(make_prediction = TRUE, tree_control =
                                        TREE.control(tree_type  = tree_type, 
                                                     d = tree.depths[dd], 
                                                     num_dir = num_dir), 
                                      shrinkage = shrinkage, nknot = nknot, loss = loss))
  val_error_vec[dd] <-  
    model_TFBoost_list[[dd]]$err_val[model_TFBoost_list[[dd]]$early_stop]
}
model_TFBoost <-  model_TFBoost_list[[which.min(val_error_vec)]]
```

We fit `fgam` which estimates the regression function using
tensor-product B-splines with roughness penalties. We set the marginal
basis to be cubic B-splines of dimension 15 and used REML to select the
parameter that penalizes the second order marginal differences.

``` r
library(refund)
nbasis_FGAM <- 15
xx <- xtrain
yy <- ytrain
model_FGAM <- fgam(yy ~ af(xx, splinepars = list(k=c(nbasis_FGAM,nbasis_FGAM), 
                                                 m = list(c(2,2),c(2,2)))), 
                   method ="REML")
```

For `TFBoost`, the prediction at early stopping time using `x_test` as
the predictor is stored in the `f_test_t` entry of the returned objects.
The returned `err_test` stored a vector of mean-squared-errors on the
test set before and at the early stopping iteration. A sanity check:

``` r
mse_TFBoost.1 <- mean((model_TFBoost$f_test_t - ytest)^2)
mse_TFBoost.2 <- model_TFBoost$err_test[model_TFBoost$early_stop]
all.equal(mse_TFBoost.1, mse_TFBoost.2)
```

    ## [1] TRUE

We compare the test errors of `TFBoost` and `fgam` and note that
`TFBoost` yields better predictions on the test set:

``` r
mse_FGAM <-  mean((predict(model_FGAM,newdata=list(xx =xtest),type='response')- ytest)^2)
print(c(mse_TFBoost.1,mse_FGAM))
```

    ## [1] 3524225 5026891

We observe that electricity demand shows a seasonality pattern likely
due to higher usage in summer for cooling and in winter for heating. To
adjust for the seasonal effects, we considered adding a
`day_of_the_year` variable as a scalar predictor, which is defined as
the number of days from January 1st of the year in which the data was
collected.

We compute the `day_of_the_year` variable and construct its values for
the training, validation, and test sets.

``` r
Year <- substr(GED$date,1,4)
day_of_the_year <- as.numeric(GED$date- as.Date(paste0("01.01.", Year), "%d.%m.%Y"))
ztrain <- matrix(day_of_the_year[idx_train], dimnames = list(NULL, "z"))
zval <-  matrix(day_of_the_year[idx_val],dimnames = list(NULL, "z"))
ztest <-  matrix(day_of_the_year[idx_test],dimnames = list(NULL, "z"))
```

Below we train `TFBoost` and `fgam` with `price` and `day_of_the_year`
as predictors:

``` r
tree.depths <- 1:4
model_TFBoost_day_list <- vector('list', length(tree.depths))
val_error_day_vec <- rep(NA, length(tree.depths))
for(dd in 1:length(tree.depths)){
  model_TFBoost_day_list[[dd]] <- 
    TFBoost(x_train = xtrain, z_train = ztrain, y_train = ytrain, x_val = xval, 
            z_val = zval,  y_val = yval, x_test = xtest, z_test = ztest, y_test = ytest,
            grid = gg, t_range  = tt, niter = niter, 
            control = TFBoost.control(make_prediction = TRUE, 
                                      tree_control = TREE.control(tree_type  = tree_type, 
                                                                  d = tree.depths[dd],
                                                                  num_dir = num_dir), 
                                      shrinkage = shrinkage, nknot = nknot, loss = loss))
  val_error_day_vec[dd] <-
    model_TFBoost_day_list[[dd]]$err_val[model_TFBoost_day_list[[dd]]$early_stop]
}
model_TFBoost_day <-  model_TFBoost_day_list[[which.min(val_error_day_vec)]]

zz <- ztrain
model_FGAM_day <- fgam(yy ~ af(xx, splinepars = list(k=c(nbasis_FGAM,nbasis_FGAM),
                                                     m = list(c(2,2),c(2,2)))) + 
                         s(zz, bs = "cs"), method = "REML")
```

We compare their test errors:

``` r
mse_TFBoost_day <-mean( (model_TFBoost_day$f_test_t- ytest)^2 )
mse_FGAM_day <-  mean( (predict(model_FGAM_day,
                               newdata=list(xx =xtest, zz = ztest),
                               type='response') - ytest)^2 )
print(c(mse_TFBoost_day,mse_FGAM_day))
```

    ## [1] 2872263 4209165

We note that the predictive performance of `TFBoost` and `fgam` both
improved comapred to only using `price` as the predictor, and `TFBoost`
yields better test errors.

For `TFBoost`, we can also separate the process of training the
predictor and evaluating it on a test set. In this way we have the
flexibility to compute predictions on multiple test sets when needed. In
this case the user needs to use `save_tree = TRUE` in the argument
`control`. This option includes the tree objects at all iterations in
the returned object, which are needed for calculating predictions. The
following code illustrates this procedure:

``` r
model_TFBoost_separate <- TFBoost(x_train = xtrain,  y_train = ytrain, 
                                  x_val = xval,  y_val = yval, 
                                  grid = gg, t_range  = tt, niter = niter, 
                                  control = 
                                    TFBoost.control(make_prediction = FALSE,
                                                    tree_control = TREE.control(
                                                      tree_type  = tree_type, 
                                                      d = which.min(val_error_vec), 
                                                      num_dir = num_dir),
                                                    shrinkage = shrinkage, 
                                                    nknot = nknot, loss = loss, 
                                                    save_tree = TRUE))
predictions <- TFBoost.predict(model_TFBoost_separate, newx = xtest)
model_TFBoost_separate <- TFBoost(x_train = xtrain, z_train = ztrain, y_train = ytrain, 
                                  x_val = xval, z_val = zval,  y_val = yval, 
                                  grid = gg, t_range  = tt, niter = niter, 
                                  control = 
                                    TFBoost.control(make_prediction = FALSE, 
                                                    tree_control = TREE.control(
                                                      tree_type  = tree_type, 
                                                      d = which.min(val_error_day_vec),
                                                      num_dir = num_dir),
                                                    shrinkage = shrinkage, 
                                                    nknot = nknot, loss = loss, 
                                                    save_tree = TRUE))
predictions_day <- TFBoost.predict(model_TFBoost_separate, newx = xtest, newz = ztest)
```

A sanity check

``` r
all.equal(mean( (predictions - ytest)^2), mse_TFBoost.1)
```

    ## [1] TRUE

``` r
all.equal(mean( (predictions_day - ytest)^2), mse_TFBoost_day)
```

    ## [1] TRUE
