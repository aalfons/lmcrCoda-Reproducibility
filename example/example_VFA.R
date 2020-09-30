# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


# load functions and packages
source("code/lmcrCoda.R")
library("VIM")
library("microbenchmark")


## check R version
if (getRversion() != "3.6.3") {
  warning("Results in the paper were obtained using R version 3.6.3.\n",
          "Your results may differ due the use of a different R version.",
          call. = FALSE)
}

## version and seed of random number generator to be used
RNGversion("3.6.3")
seed <- 20190930


# load data
load("example/VFA.RData")
dim(VFA)

# extract compositional data
X <- VFA[, 4:9]
n <- nrow(X)
D <- ncol(X)
apply(is.na(X), 2, any)  # check if there are any missing values
apply(X, 2, min)         # check that there are no zeros

# scale compositional observations
X <- 1000 * X / rowSums(X)

# real-valued covariates
R <- log(VFA[, 10:12])
p <- ncol(R)
names(R) <- paste0("ln", names(R))

# dummy variable
Diet <- VFA[, 2]
Dummies <- data.frame(Diet = as.integer(Diet == "Mixed"))

# response variable
y <- log(VFA[, 3])

# control parameters
tau <- 0.99
numDiscrete <- 10
corrlim <- NA
control <- lmrob.control(max.it = 1000, k.max = 10000, maxit.scale = 1000)

## our method
set.seed(seed)
fitBFMI <- lmcrCoda(X, y, R, Dummies, tau = tau, numDiscrete = numDiscrete,
                    control = control, corrlim = corrlim)

# fit OLS and MM with CoDa
covariates <- cbind(R, DietMixed = as.numeric(Diet == "Mixed"))
fitOLS <- lmCoda(X, y, covariates = covariates)
fitMM <- lmrobCoda(X, y, covariates = covariates, control = control)

# compare methods
fitOLS
fitMM
fitBFMI


# prepare data set with filtered cells for plotting
VFAfiltered <- cbind(logCH4 = setNA(y, fitBFMI$indicesY),
                     mapply(setNA, X, fitBFMI$indicesX),
                     mapply(setNA, R, fitBFMI$indicesR))
VFAfiltered[fitBFMI$indicesRow, ] <- NA
# nicer labels for plot
labels <- c(sprintf("ln(%s)", names(VFA)[3]), names(VFA)[4:9],
            sprintf("ln(%s)", names(VFA)[10:12]))
# create plot of filtered cells
pdf("example/VFA_filtered.pdf", width = 6, height = 6)
par(las = 1, mar = c(5.5, 4, 0.5, 0) + 0.1)
matrixplot(VFAfiltered, labels = labels)
dev.off()


# compute percentage of outliers
nOutRow <- length(fitBFMI$indicesRow)
nOutCell <- sum(sapply(fitBFMI$indicesX, length)) +
  sum(sapply(fitBFMI$indicesR, length)) + length(fitBFMI$indicesY)
100 * nOutRow / n
100 * nOutCell / ((n-nOutRow)*(D+p+2))


# sensitivity analysis
taus <- seq(from = 0.995, to = 0.975, by = -0.005)
# without variable selection in imputations
fitList1 <- lapply(taus, function(tau) {
  set.seed(seed)
  lmcrCoda(X, y, R, Dummies, tau = tau, numDiscrete = numDiscrete,
                   control = control, corrlim = NA)
})
# with variable selection in imputations
fitList2 <- lapply(taus, function(tau) {
  set.seed(seed)
  lmcrCoda(X, y, R, Dummies, tau = tau, numDiscrete = numDiscrete,
           control = control, corrlim = 0.2)
})
# add names giving values of tau and show results
names(fitList1) <- names(fitList2) <- as.character(taus)
fitList1
fitList2


# check computation time
microbenchmark(
  OLS = lmCoda(X, y, covariates = covariates),
  MM = lmrobCoda(X, y, covariates = covariates, control = control),
  BFMI = lmcrCoda(X, y, R, Dummies, tau = tau, numDiscrete = numDiscrete,
                  control = control, corrlim = corrlim),
  times = 10
)
