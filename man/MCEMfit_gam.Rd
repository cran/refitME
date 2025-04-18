% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refitME_package.r
\name{MCEMfit_gam}
\alias{MCEMfit_gam}
\title{Function for wrapping the MCEM algorithm on \code{gam} objects}
\source{
See \url{https://github.com/JakubStats/refitME} for an RMarkdown vignette with examples. With permission from Matt Wand, we have now made these data available in the refitME R-package.
}
\usage{
MCEMfit_gam(
  mod,
  family,
  sigma.sq.u,
  B = 50,
  epsilon = 1e-05,
  silent = FALSE,
  ...
)
}
\arguments{
\item{mod}{: a \code{gam} object (this is the naive fitted model). Make sure the first \eqn{p} input predictor variables entered in the naive model are the specified error-contaminated variables. These \eqn{p} predictors also need the measurement error variance to be specified in \code{sigma.sq.u}, see below.}

\item{family}{: a specified family/distribution.}

\item{sigma.sq.u}{: measurement error (ME) variance. A scalar if there is only one error-contaminated predictor variable, otherwise this must be stored as a vector (of ME variances) or a matrix if the ME covariance matrix is known.}

\item{B}{: the number of Monte Carlo replication values (default is set to 50).}

\item{epsilon}{: convergence threshold (default is set to 0.00001).}

\item{silent}{: if \code{TRUE}, the convergence message (which tells the user if the model has converged and reports the number of iterations required) is suppressed (default is set to \code{FALSE}).}

\item{...}{: further arguments passed to \code{gam}.}
}
\value{
\code{MCEMfit_gam} returns the original naive fitted model object but coefficient estimates and the covariance matrix have been replaced with the final MCEM model fit. Standard errors and the effective sample size (which diagnose how closely the proposal distribution matches the posterior, see equation (2) of Stoklosa, Hwang and Warton) have also been included as outputs.
}
\description{
Function for wrapping the MCEM algorithm on GAMs where predictors are subject to measurement error/error-in-variables.
}
\examples{
# A GAM example using the air pollution data set from the SemiPar package.

library(refitME)
library(mgcv)
library(dplyr)

data(Milanmortdata)

dat.air <- sample_n(Milanmortdata, 100) # Takes a random sample of size 100.

Y <- dat.air[, 6]  # Mortality counts.

n <- length(Y)

z1 <- (dat.air[, 1])
z2 <- (dat.air[, 4])
z3 <- (dat.air[, 5])
w1 <- log(dat.air[, 9])  # The error-contaminated predictor (total suspended particles).

dat <- data.frame(cbind(Y, w1, z1, z2, z3))

gam_naiv <- gam(Y ~ s(w1), family = "poisson", data = dat)

sigma.sq.u <- 0.0915 # Measurement error variance.

B <- 10  # Consider increasing this if you want a more accurate answer.

gam_MCEM <- refitME(gam_naiv, sigma.sq.u, B)

plot(gam_MCEM, select = 1)

detach(package:mgcv)

}
\references{
Ganguli, B, Staudenmayer, J., and Wand, M. P. (2005). Additive models with predictors subject to measurement error. \emph{Australian & New Zealand Journal of Statistics}, \strong{47}, 193–202.

Stoklosa, J., Hwang, W-H., and Warton, D.I. \pkg{refitME}: Measurement Error Modelling using Monte Carlo Expectation Maximization in \proglang{R}.
}
\seealso{
\code{\link{MCEMfit_glm}}
}
\author{
Jakub Stoklosa, Wen-Han Hwang and David I. Warton.
}
