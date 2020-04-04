# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


## load required packages
library("mvtnorm")
library("robCompositions")

## function to compute pivot coordinates based on available information
pivotCoordNA <- function(x, pivotvar = 1) {
  # initializations
  n <- nrow(x)
  D <- ncol(x)
  # order according to pivot variable
  order <- c(pivotvar, setdiff(seq_len(D), pivotvar))
  x <- x[, order]
  # transform observations
  z_list <- lapply(seq_len(n), function(i) {
    # find observed parts
    xi <- x[i, ]
    observed <- !is.na(xi)
    D_obs <- sum(observed)
    # initialize the pivot coordinates
    zi <- rep.int(NA_real_, D-1)
    j <- 0
    # loop over compositional parts and update scaling factor according to
    # available information
    for (k in seq_len(D-1)) {
      if (observed[k] && any(observed[(k+1):D])) {
        j <- j + 1
        zi[k] <- - sqrt((D_obs-j)/(D_obs-j+1)) * log(gm(xi[(k+1):D]) / xi[k])
      }
    }
    if (j != D_obs-1) {
      warning(sprintf("something want wrong: D_obs = %d, j = %d", D_obs, j))
    }
    # return pivot coordinates
    zi
  })
  # combine observations into matrix
  z <- do.call(rbind, z_list)
  # add column names
  if(all(nchar(colnames(x)) > 1)) {
    colnames(z) <- sapply(seq_len(D-1), function(j, labels) {
      denominator <- paste(substr(labels[(j+1):D], 1, 2), collapse = "-")
      paste(labels[j], "_", denominator, sep = "")
    }, labels = colnames(x))
  }
  # return pivot coordinates
  z
}

## function for linear regression based on covariance matrix with pairwise
## complete observations (to handle missing values directly on the simplex)
lmNA <- function(x, y) {
  # initializations
  z <- cbind(y, x)
  n <- nrow(z)
  # compute covariance matrix
  # (an eigenvalue correction if necessary, but there were no issues here)
  S <- cov(z, use = "pairwise.complete.obs")
  # compute regression coefficients
  b <- drop(solve(S[-1, -1]) %*% S[-1, 1])
  a <- mean(y, na.rm = TRUE) - sum(colMeans(x, na.rm = TRUE) * b)
  # return coefficients
  c("(Intercept)" = a, b)
}


## check R version
if (getRversion() != "3.6.3") {
  warning("Results in the paper were obtained using R version 3.6.3.\n",
          "Your results may differ due the use of a different R version.",
          call. = FALSE)
}

## version and seed of random number generator to be used
RNGversion("3.6.3")
seed <- 20190702


## control parameters
n <- 250    # number of observations
D <- 6      # number of compositional parts
pNA <- 0.1  # probability of having a missing value
# covariance matrix of the explanatory variables on the simplex
Sigma <- 0.5^t(sapply(1:(D-1), function(i, j) abs(i-j), 1:(D-1)))
# regression coefficients on the simplex
beta <- rep(c(1, 0), length.out = D-1)
# number of simulation runs
R <- 500


## run simulation
set.seed(seed)
cat(paste(Sys.time(), ": starting ...\n"))
results <- lapply(1:R, function(r) {

  ## print simulation run
  cat(paste(Sys.time(), sprintf(":   run = %d\n", r)))

  ## generate data
  # complete data
  z <- rmvnorm(n, sigma = Sigma)
  y <- drop(z %*% beta) + rnorm(n)
  x <- pivotCoordInv(z)
  colnames(x) <- paste0("X", 1:D)
  # replace some missing values (MCAR)
  xNA <- x
  xNA[runif(n*D) < pNA] <- NA_real_

  # ensure that an error doesn't stop the entire simulation
  # (in case of an error, NULL is returned)
  tryCatch({

    ## compute regression coefficients for different pivot coordinates
    # first variable as first coordinate
    z1 <- pivotCoordNA(xNA, pivotvar = 1)
    coef1 <- lmNA(z1, y)
    # second variable as first coordinate
    z2 <- pivotCoordNA(xNA, pivotvar = 2)
    coef2 <- lmNA(z2, y)

    ## compare out of sample predictions, which should be identical for
    ## different pivot coordinates
    # we only need test data for the pivot coordinates (we don't need
    # to generate the response, we're only interested in comparing the
    # predictions based on two sets of pivot coordinates with each other)
    zTest <- rmvnorm(n, sigma = Sigma)
    xTest <- pivotCoordInv(zTest)
    # generate pivot coordinates
    zTest1 <- pivotCoord(xTest, pivotvar = 1)
    zTest2 <- pivotCoord(xTest, pivotvar = 2)
    # compute predictions
    yHat1 <- coef1[1] + drop(as.matrix(zTest1) %*% coef1[-1])
    yHat2 <- coef2[1] + drop(as.matrix(zTest2) %*% coef2[-1])
    # return average of absolute difference between predictions
    mean(abs(yHat1 - yHat2))

  }, error = function(e) NULL, warning = function(w) NULL)

})

## combine results into vector
results <- unlist(results)
cat(paste(Sys.time(), ": finished.\n"))

# summarize results
m <- mean(results)
s <- sd(results)
cat(sprintf("\nMean = %.4f, Std. Dev. = %.4f\n", m, s))

## store results
file <- "simulations/results/results_pivotCoordNA.RData"
save(results, file = file)
