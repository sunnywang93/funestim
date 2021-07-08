################################################################################
#                   Functions for covariance estimation                        #
################################################################################

#' Perform an estimation of the covariance with smoothing splines.
#' 
#' This function performs the estimation of the covariance of a set of curves
#' using smoothing splines.
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  }
#' @param U Vector, sampling points at which estimate the covariance.
#' @param nbasis Integer (default=5), number of basis to use for the splines.
#' @param centered Boolean (default=FALSE), are the data centered?
#' @param nodiag Boolean (default=TRUE), should the diagonal be removed from
#'  the fitting of the model?
#' 
#' @return A matrix representing the covariance surface.
#' 
#' @references Cai, T., Yuan, M. (2010) - Nonparametric covariance function 
#' estimation for functional and longitudinal data. University of Pennsylvania 
#' and Georgia inistitute of technology
#' @export
covariance_ss <- function(data, U, nbasis = 5, centered = FALSE, nodiag = TRUE) {
  predict_ssanova <- utils::getFromNamespace("predict.ssanova", "gss")
  
  if(!inherits(data, 'list')) data <- checkData(data)
  data_ <- list2cai(data)
  time <- data_$time
  x <- data_$x
  subject <- data_$obs
  
  if (!centered) {
    fit <- stats::smooth.spline(time, x)
    x <- x - stats::fitted(fit)
  }
  gg <- NULL
  for(zz in unique(subject)) {
    if(sum(subject == zz) > 1) {
      tt <- time[subject == zz]
      xx <- x[subject == zz]
      g <- expand.grid(t1 = tt, t2 = tt)
      scov <- xx %*% t(xx)
      if(nodiag) scov <- scov + diag(rep(Inf, length(xx)))
      g$z <- matrix(scov, ncol = 1)
      gg <- rbind(gg, g[g$z < Inf, ])
    }
  }
  
  gg <- unique(gg)
  nobs <- nrow(gg)
  tt <- min(time) + (max(time) - min(time)) * (1:nbasis)/(nbasis + 1)
  g <- expand.grid(t1 = tt, t2 = tt)
  g$z <- 0
  gg <- rbind(g, gg)
  
  fit <- gss::ssanova(z ~ t1 * t2, data = gg, id.basis = 1:(nbasis * nbasis))
  
  new <- expand.grid(t1 = U, t2 = U)
  estim <- predict_ssanova(fit, newdata = new)
  matrix(estim, ncol = length(U), nrow = length(U))
}

#' Perform an estimation of the covariance with local linear smoothers.
#' 
#' This function performs the estimation of the covariance of a set of curves 
#' using local linear smoothers where the bandwidth is given.
#' 
#' @param data A list, where each element represents a curve. Each curve have to
#'  be defined as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points.
#'  } 
#' @param U Vector, sampling points at which estimate the covariance.
#' @param b Numeric (default=0.1), the bandwidth to use.
#' 
#' @return A matrix representing the covariance surface.
#' 
#' @references Zhang X. and Wang J.-L. (2016) - From sparse to dense functional
#'  data and beyond, The Annals of Statistics
#' @export
covariance_lll <- function(data, U, b = 0.1) {
  if(!inherits(data, 'list')) data <- checkData(data)
  data_ <- list2cai(data)
  L3 <- fdapace::MakeFPCAInputs(IDs = data_$obs, 
                                tVec = data_$time, 
                                yVec = data_$x)
  fdapace::GetCovSurface(L3$Ly, L3$Lt, 
                         list(kernel = 'epan',
                              nRegGrid = length(U),
                              methodMuCovEst = 'smooth',
                              userBwCov = b))$cov
}
