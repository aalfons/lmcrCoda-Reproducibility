# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


## The following packages are required by lmcrCoda() and shootingS().  They
## are included here as well to make packrat work properly.
library("cellWise")
library("robCompositions")
library("robustHD")

## load required functions and packages
devtools::source_url(
  "https://raw.github.com/aalfons/lmcrCoda/master/lmcrCoda.R"
)
devtools::source_url(
  "https://raw.github.com/aalfons/shootingS/master/functions_shootingS.R"
)
library("microbenchmark")
library("mvtnorm")
library("robreg3S")
options(warn = 1)  # show warnigns immediately


## check R version
if (getRversion() != "3.6.3") {
  warning("Results in the paper were obtained using R version 3.6.3.\n",
          "Your results may differ due the use of a different R version.",
          call. = FALSE)
}

## version and seed of random number generator to be used
RNGversion("3.6.3")
seed <- 20190813


## control parameters for simulation
dims <- expand.grid(n = c(50, 100, 200), D = c(5, 10, 20))
R <- 100          # number of simulation runs

## control parameters for contamination
zeta <- 0.02      # contamination level
shift <- 5        # shift for leverage points
multiplier <- 10  # multiplier for cellwise outliers

## control parameters for methods
tau <- 0.99
pOutLR <- 0.5
pOutRow <- 0.75
lmrobControl <- lmrob.control(k.max = 100000, max.it = 5000,
                              maxit.scale = 5000)
corrlim <- NULL


## run simulation
cat(paste(Sys.time(), ": starting ...\n"))
results_list <- mapply(function(n, D) {

  ## use the same seed for each dimension
  ## (then simulations can be run separately for different dimensions)
  set.seed(seed)

  ## print dimensions
  cat(paste(Sys.time(), sprintf(": n = %d, D = %d\n", n, D)))

  ## further parameters for data generation
  # covariance matrix of the explanatory variables on the simplex
  Sigma0 <- 0.5^t(sapply(1:(D-1), function(i, j) abs(i-j), 1:(D-1))) / 10
  loadings <- eigen(Sigma0)$vectors
  # regression coefficients (without intercept)
  beta <- rep_len(1:0, length.out = D-1)  # for good data points
  gamma <- rep.int(-1, times = D-1)       # for outliers
  # regression coefficients (with intercept)
  alphabeta <- c(0, beta)
  # standard deviation of error terms
  sigma <- 0.25

  ## perform simulation for current sample size
  results_nD <- lapply(1:R, function(r) {

    ## print simulation run
    cat(paste(Sys.time(), sprintf(":   run = %d\n", r)))

    ## generate pivot coordinates and error terms
    Z <- rmvnorm(n, sigma = Sigma0)
    e <- rnorm(n, sd = sigma)

    ## indices of cellwise and rowwise outliers
    # first determine outlying rows
    indicesRow <- (runif(n) < zeta)
    nOK <- n - sum(indicesRow)
    # determine outlying cells in non-outlying rows of predictor matrix
    indicesX <- matrix(FALSE, nrow = n, ncol = D)
    indicesX[!indicesRow] <- (runif(nOK * D) < zeta)
    # determine outlying cells in non-outlying rows of response
    indicesY <- rep.int(FALSE, n)
    indicesY[!indicesRow] <- (runif(nOK) < zeta)

    ## generate response
    y <- drop(Z %*% beta) + e

    ## generate rowwise outliers
    # put leverage points in smallest principal component
    Q <- Z %*% loadings
    Q[indicesRow, D-1] <- Q[indicesRow, D-1] + shift
    Z[indicesRow, ] <- Q[indicesRow, ] %*% t(loadings)
    # genreate response with very different regression coefficients
    y[indicesRow] <- drop(Z[indicesRow, ] %*% gamma) + e[indicesRow]

    ## transform to compositional space
    X <- pivotCoordInv(Z)

    ## generate cellwise outliers
    X[indicesX] <- multiplier * X[indicesX]
    y[indicesY] <- multiplier * y[indicesY]

    ## measure computation time of different methods
    ## (Transforming the compositions to the required coordinate system is
    ## considered part of the algorithm and therefore included in measuring
    ## the computation time.  This is done within the function for some
    ## methods, but other functions are more general and require the
    ## compositions to be transformed beforehand.)

    # compositional OLS
    dfOLS <- tryCatch({
      cpu <- microbenchmark(OLS = lmCoda(X, y), times = 1)
      df_cpu <- summary(cpu, unit = "s")[, c("expr", "mean")]
      names(df_cpu) <- c("Method", "CPU")
      data.frame(n = n, D = D, Run = r, df_cpu)
    }, error = function(e) NULL, warning = function(w) NULL)
    # compositional MM-estimator
    dfMM <- tryCatch({
      cpu <- microbenchmark(MM = lmrobCoda(X, y, control = lmrobControl),
                            times = 1)
      df_cpu <- summary(cpu, unit = "s")[, c("expr", "mean")]
      names(df_cpu) <- c("Method", "CPU")
      data.frame(n = n, D = D, Run = r, df_cpu)
    }, error = function(e) NULL, warning = function(w) NULL)
    # proposed method
    dfBF <- tryCatch({
      cpu <- microbenchmark(
        "BF-MI" = lmcrCoda(X, y, tau = tau, pOutLR = pOutLR,
                           pOutRow = pOutRow, control = lmrobControl,
                           corrlim = corrlim, interpretable = TRUE,
                           single = FALSE),
        times = 1)
      df_cpu <- summary(cpu, unit = "s")[, c("expr", "mean")]
      names(df_cpu) <- c("Method", "CPU")
      data.frame(n = n, D = D, Run = r, df_cpu)
    }, error = function(e) NULL, warning = function(w) NULL)
    # shooting S-estimator with Tukey biweight loss function
    dfShSBi <- tryCatch({
      cpu <- microbenchmark(
        "shS-bi" = {
          D <- ncol(X)
          combinations <- combn(seq_len(D), 2, simplify = FALSE)
          Xlr <- lapply(combinations, function(j) log(X[, j[1]] / X[, j[2]]))
          Xlr <- do.call(cbind, Xlr)
          fitShSBi <- shooting(Xlr, y, method = "biweight")
        }, times = 1)
      df_cpu <- summary(cpu, unit = "s")[, c("expr", "mean")]
      names(df_cpu) <- c("Method", "CPU")
      data.frame(n = n, D = D, Run = r, df_cpu)
    }, error = function(e) NULL, warning = function(w) NULL)
    # 3-step regression
    df3S <- tryCatch({
      cpu <- microbenchmark(
        "3S" = {
          D <- ncol(X)
          denominator <- sample.int(D, 1)
          Xalr <- log(X[, -denominator] / X[, denominator])
          robreg3S(y, Xalr, maxiter = 5000)
        }, times = 1)
      df_cpu <- summary(cpu, unit = "s")[, c("expr", "mean")]
      names(df_cpu) <- c("Method", "CPU")
      data.frame(n = n, D = D, Run = r, df_cpu)
    }, error = function(e) NULL, warning = function(w) NULL)

    # results for current simulation run
    rbind(dfOLS, dfMM, dfBF, dfShSBi, df3S)

  })

  ## combine results for current dimensions into data frame
  do.call(rbind, results_nD)

}, n = dims$n, D = dims$D, SIMPLIFY = FALSE, USE.NAMES = FALSE)

## combine results into data frame
results <- do.call(rbind, results_list)
cat(paste(Sys.time(), ": finished.\n"))

## store results
file <- "simulations/results/results_cpu_time.RData"
save(results, seed, file = file)
