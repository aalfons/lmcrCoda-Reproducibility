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
library("mvtnorm")
library("robreg3S")
options(warn = 1)  # show warnigns immediately

## additional functions
MSE <- function(true, estimate) mean((true - estimate)^2)


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
ks <- 1:3                              # multiplication of variation matrix
R <- 1000                              # number of simulation runs

## control parameters for contamination
zetas <- c(0, 0.02, 0.05, 0.10, 0.20)  # contamination levels
shift <- 5                             # shift for leverage points
multiplier <- 10                       # multiplier for cellwise outliers

## control parameters for methods
tau <- 0.99
pOutLR <- 0.5
pOutRow <- 0.75
lmrobControl <- lmrob.control(k.max = 100000, max.it = 5000,
                              maxit.scale = 5000)
corrlim <- NULL


## run simulation
## (assingment to object is just to prevent printing NULL,
## results are saved to disk within loop)
cat(paste(Sys.time(), ": starting ...\n"))
results_null <- mapply(function(n, D) {

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

  ## pairwise combinations of variables (for computing logratios)
  combinations <- combn(seq_len(D), 2, simplify = FALSE)

  ## perform simulation for current sample size
  results_nD <- lapply(1:R, function(r) {

    ## print simulation run
    cat(paste(Sys.time(), sprintf(":   run = %d\n", r)))

    ## generate pivot coordinates and error terms (use the same as basis for
    ## all multipliers of the variation matrix for maximum comparability)
    # data for model fitting
    Z <- rmvnorm(n, sigma = Sigma0)
    e <- rnorm(n, sd = sigma)
    # independent test data for prediction
    ZTest <- rmvnorm(n, sigma = Sigma0)
    eTest <- rnorm(n, sd = sigma)

    ## loop over contamination levels
    results_nDr <- lapply(zetas, function(zeta) {

      # ## print contamination level
      # cat(paste(Sys.time(), sprintf(":     zeta = %.2f\n", zeta)))

      ## indices of cellwise and rowwise outliers
      ## Ideally, we would draw from a uniform distribution to determine which
      ## rows/cells are outliers and use the same random numbers for all
      ## contamination levels for maximum comparability.  However, that doesn't
      ## work if we mix rowwise and cellwise outliers such that cellwise outliers
      ## are generated only in non-outlying rows.  But we can still use the same
      ## rows/cells across multipliers of the variation matrix.
      # first determine outlying rows
      indicesRow <- (runif(n) < zeta)
      nOK <- n - sum(indicesRow)
      # determine outlying cells in non-outlying rows of predictor matrix
      indicesX <- matrix(FALSE, nrow = n, ncol = D)
      indicesX[!indicesRow] <- (runif(nOK * D) < zeta)
      # determine outlying cells in non-outlying rows of response
      indicesY <- rep.int(FALSE, n)
      indicesY[!indicesRow] <- (runif(nOK) < zeta)

      ## loop over different multipliers of the variation matrix
      results_nDrc <- lapply(ks, function(k) {

        # ## print multiplier
        # cat(paste(Sys.time(), sprintf(":       k = %d\n", k)))

        ## generate data
        # multiply the pivot coordinates with the square root of the multiplier
        Z <- sqrt(k) * Z
        # generate response
        y <- drop(Z %*% beta) + e

        ## generate rowwise outliers
        # put leverage points in smallest principal component
        # (make sure that the shift scales with the variation matrix)
        Q <- Z %*% loadings
        Q[indicesRow, D-1] <- Q[indicesRow, D-1] + shift * sqrt(k)
        Z[indicesRow, ] <- Q[indicesRow, ] %*% t(loadings)
        # genreate response with very different regression coefficients
        y[indicesRow] <- drop(Z[indicesRow, ] %*% gamma) + e[indicesRow]

        ## transform to compositional space
        X <- pivotCoordInv(Z)

        ## generate cellwise outliers
        X[indicesX] <- multiplier * X[indicesX]
        y[indicesY] <- multiplier * y[indicesY]

        ## obtain transformations
        # compute pivot coordinates
        Z <- as.matrix(pivotCoord(X))
        # compute pairwise logratios
        Xlr <- lapply(combinations, function(j) log(X[, j[1]] / X[, j[2]]))
        Xlr <- do.call(cbind, Xlr)
        # compute alr coordinates with randomly selected denominator
        denominator <- sample.int(D, 1)
        Xalr <- log(X[, -denominator] / X[, denominator])

        ## generate test data for prediction evaluation
        # multiply the pivot coordinates with the square root of the multiplier
        ZTest <- sqrt(k) * ZTest
        # generate response
        yTest <- drop(ZTest %*% beta) + eTest
        # transform to compositional space
        XTest <- pivotCoordInv(ZTest)
        # compute pairwise logratios
        XlrTest <- lapply(combinations, function(j) {
          log(XTest[, j[1]] / XTest[, j[2]])
        })
        XlrTest <- do.call(cbind, XlrTest)
        # compute alr coordinates
        XalrTest <- log(XTest[, -denominator] / XTest[, denominator])

        ## OLS and MM
        predictors <- cbind("(Intercept)" = rep.int(1, n), Z)
        predictorsTest <- cbind("(Intercept)" = rep.int(1, n), ZTest)
        # apply OLS estimator
        dfOLS <- tryCatch({
          fitOLS <- lm.fit(predictors, y)              # fit model
          coefOLS <- coef(fitOLS)                      # extract coefficients
          yHatOLS <- drop(predictorsTest %*% coefOLS)  # compute predictions
          data.frame(n = n, D = D, Run = r, Zeta = zeta, k = k,
                     Method = "OLS", MSE = MSE(alphabeta, coefOLS),
                     MSEP = MSE(yTest, yHatOLS), stringsAsFactors = FALSE)
        }, error = function(e) NULL, warning = function(w) NULL)
        # apply MM-estimator
        dfMM <- tryCatch({
          fitMM <- lmrob.fit(predictors, y, control = lmrobControl)
          coefMM <- coef(fitMM)                      # extract coefficients
          yHatMM <- drop(predictorsTest %*% coefMM)  # compute predictions
          data.frame(n = n, D = D, Run = r, Zeta = zeta, k = k,
                     Method = "MM", MSE = MSE(alphabeta, coefMM),
                     MSEP = MSE(yTest, yHatMM), stringsAsFactors = FALSE)
        }, error = function(e) NULL, warning = function(w) NULL)

        ## proposed method with ideal filter of true outliers
        dfIF <- tryCatch({
          # apply method
          fitIF <- lmcrCoda(X, y, tau = tau, pOutLR = pOutLR,
                            pOutRow = pOutRow, control = lmrobControl,
                            corrlim = corrlim, interpretable = FALSE,
                            single = TRUE, indicesX = indicesX,
                            indicesY = indicesY, indicesRow = indicesRow)
          # extract coefficients
          coefIFSI <- fitIF$SI[, 1]
          coefIFMI <- fitIF$MI[, 1]
          # compute MSE
          mseIFSI <- MSE(alphabeta, coefIFSI)
          mseIFMI <- MSE(alphabeta, coefIFMI)
          # compute predictions
          yHatIFSI <- drop(predictorsTest %*% coefIFSI)
          yHatIFMI <- drop(predictorsTest %*% coefIFMI)
          # compute MSEP
          msepIFSI <- MSE(yTest, yHatIFSI)
          msepIFMI <- MSE(yTest, yHatIFMI)
          # combine into data frame
          data.frame(n = n, D = D, Run = r, Zeta = zeta,
                     k = k, Method = c("IF-SI", "IF-MI"),
                     MSE = c(mseIFSI, mseIFMI),
                     MSEP = c(msepIFSI, msepIFMI),
                     stringsAsFactors = FALSE)
        }, error = function(e) NULL, warning = function(w) NULL)

        ## proposed method
        dfBF <- tryCatch({
          fitBF <- lmcrCoda(X, y, tau = tau, pOutLR = pOutLR,
                            pOutRow = pOutRow, control = lmrobControl,
                            corrlim = corrlim, interpretable = FALSE,
                            single = TRUE)
          # extract coefficients
          coefBFSI <- fitBF$SI[, 1]
          coefBFMI <- fitBF$MI[, 1]
          # compute MSE
          mseBFSI <- MSE(alphabeta, coefBFSI)
          mseBFMI <- MSE(alphabeta, coefBFMI)
          # compute predictions
          yHatBFSI <- drop(predictorsTest %*% coefBFSI)
          yHatBFMI <- drop(predictorsTest %*% coefBFMI)
          # compute MSEP
          msepBFSI <- MSE(yTest, yHatBFSI)
          msepBFMI <- MSE(yTest, yHatBFSI)
          # combine into data frame
          data.frame(n = n, D = D, Run = r, Zeta = zeta,
                     k = k, Method = c("BF-SI", "BF-MI"),
                     MSE = c(mseBFSI, mseBFMI),
                     MSEP = c(msepBFSI, msepBFMI),
                     stringsAsFactors = FALSE)
        }, error = function(e) NULL, warning = function(w) NULL)

        ## shooting S-estimator on pairwise logratios
        predictorsTest <- cbind("(Intercept)" = rep.int(1, n), XlrTest)
        # Tukey biweight loss function
        dfShSBi <- tryCatch({
          fitShSBi <- shooting(Xlr, y, method = "biweight")
          coefShSBi <- fitShSBi$coef                      # extract coefficients
          yHatShSBi <- drop(predictorsTest %*% coefShSBi) # compute predictions
          data.frame(n = n, D = D, Run = r, Zeta = zeta, k = k,
                     Method = "shS-bi", MSE = NA_real_,
                     MSEP = MSE(yTest, yHatShSBi),
                     stringsAsFactors = FALSE)
        }, error = function(condition) NULL)
        # skipped Huber loss function
        dfShSSkH <- tryCatch({
          fitShSSkH <- shooting(Xlr, y, method = "skHuber")
          coefShSSkH <- fitShSSkH$coef                     # extract coefficients
          yHatShSSkH <- drop(predictorsTest %*% coefShSSkH)# compute predictions
          data.frame(n = n, D = D, Run = r, Zeta = zeta, k = k,
                     Method = "shS-skH", MSE = NA_real_,
                     MSEP = MSE(yTest, yHatShSSkH),
                     stringsAsFactors = FALSE)
        }, error = function(e) NULL, warning = function(w) NULL)

        ## 3S regression on alr coordinates
        predictorsTest <- cbind("(Intercept)" = rep.int(1, n), XalrTest)
        # comppute 3S regression
        df3S <- tryCatch({
          fit3S <- robreg3S(y, Xalr, maxiter = 5000)
          coef3S <- coef(fit3S)                      # extract coefficients
          yHat3S <- drop(predictorsTest %*% coef3S)  # compute predictions
          msep3S <- MSE(yTest, yHat3S)               # compute MSEP
          data.frame(n = n, D = D, Run = r, Zeta = zeta, k = k,
                     Method = "3S", MSE = NA_real_,
                     MSEP = MSE(yTest, yHat3S),
                     stringsAsFactors = FALSE)
        }, error = function(e) NULL, warning = function(w) NULL)

        ## results for current multiplier of variation matrix
        rbind(dfOLS, dfMM, dfIF, dfBF, dfShSBi, dfShSSkH, df3S)

      })

      ## combine results for current contamination level into data frame
      do.call(rbind, results_nDrc)

    })

    ## combine results for current simulation run into data frame
    do.call(rbind, results_nDr)

  })

  ## combine results for current dimensions into data frame
  results_nD <- do.call(rbind, results_nD)

  ## store results for current sample size
  file <- "simulations/results/results_main_%d_%d.RData"
  save(results_nD, seed, file = sprintf(file, n, D))

}, n = dims$n, D = dims$D, SIMPLIFY = FALSE, USE.NAMES = FALSE)

## end simulation
cat(paste(Sys.time(), ": finished.\n"))
