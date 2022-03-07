
[![Build
Status](https://app.travis-ci.com/AnthonyChristidis/stepSplitReg.svg?branch=master)](https://app.travis-ci.com/AnthonyChristidis/stepSplitReg)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/stepSplitReg)](https://cran.r-project.org/package=stepSplitReg)
[![Downloads](http://cranlogs.r-pkg.org/badges/stepSplitReg)](https://cran.r-project.org/package=stepSplitReg)

# stepSplitReg

This package provides functions for performing stepwise split
regularized regression.

------------------------------------------------------------------------

### Installation

You can install the **stable** version on [R
CRAN](https://cran.r-project.org/package=stepSplitReg).

``` r
install.packages("stepSplitReg", dependencies = TRUE)
```

You can install the **development** version from
[GitHub](https://github.com/AnthonyChristidis/stepSplitReg)

``` r
library(devtools)
devtools::install_github("AnthonyChristidis/stepSplitReg")
```

### Usage

``` r
# Required Libraries
library(mvnfast)

# Setting the parameters
p <- 800
n <- 40
n.test <- 2000
sparsity <- 0.2
rho <- 0.5
SNR <- 3
set.seed(0)

# Generating the coefficient
p.active <- floor(p*sparsity)
a <- 4*log(n)/sqrt(n)
neg.prob <- 0.2
nonzero.betas <- (-1)^(rbinom(p.active, 1, neg.prob))*(a + abs(rnorm(p.active)))

# Correlation structure
Sigma <- matrix(0, p, p)
Sigma[1:p.active, 1:p.active] <- rho
diag(Sigma) <- 1
true.beta <- c(nonzero.betas, rep(0 , p - p.active))

# Computing the noise parameter for target SNR
sigma.epsilon <- as.numeric(sqrt((t(true.beta) %*% Sigma %*% true.beta)/SNR))

# Simulate some data
set.seed(1)
x.train <- mvnfast::rmvn(n, mu=rep(0,p), sigma=Sigma)
y.train <- 1 + x.train %*% true.beta + rnorm(n=n, mean=0, sd=sigma.epsilon)
x.test <- mvnfast::rmvn(n.test, mu=rep(0,p), sigma=Sigma)
y.test <- 1 + x.test %*% true.beta + rnorm(n.test, sd=sigma.epsilon)

# Stepwise Split Regularized Regression
step.out <- cv.stepSplitReg(x.train, y.train, n_models = c(5, 10), max_variables = NULL, keep = 4/4,
                            model_criterion = c("F-test", "RSS")[1],
                            stop_criterion = c("F-test", "pR2", "aR2", "R2", "Fixed")[1], stop_parameter = 0.05, 
                            shrinkage = TRUE, alpha = 4/4, include_intercept = TRUE, 
                            n_lambda = 50, tolerance = 1e-2, max_iter = 1e5, n_folds = 5, 
                            model_weights = c("Equal", "Proportional", "Stacking")[1], 
                            n_threads = 1)
step.coefficients <- coef(step.out, group_index = 1:20)
step.predictions <- predict(step.out, x.test, group_index = 1:20)
mspe.step <- mean((step.predictions-y.test)^2)/sigma.epsilon^2
```

### License

This package is free and open source software, licensed under GPL (&gt;=
2).
