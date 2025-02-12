% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_mean.R
\name{mean_ll}
\alias{mean_ll}
\title{Perform an estimation of the mean with local linear smoothers.}
\usage{
mean_ll(
  data,
  U = seq(0, 1, length.out = 101),
  t0_list = 0.5,
  grid = NULL,
  nb_obs_minimal = 2,
  kernel = "epanechnikov"
)
}
\arguments{
\item{data}{A list, where each element represents a curve. Each curve have to
be defined as a list with two entries:
\itemize{
\item \strong{$t} The sampling points
\item \strong{$x} The observed points.
}}

\item{U}{Vector, sampling points at which estimate the curves.}

\item{t0_list}{Vector, the sampling points at which we estimate the
parameters.}

\item{grid}{Vector, default=NULL. A grid of bandwidths.}

\item{nb_obs_minimal}{Integer, minimum number of observation for the smoothing.}

\item{kernel}{Character string, the kernel used for the estimation:
\itemize{
\item epanechnikov (default)
\item uniform
\item biweight
}}
}
\value{
A list of with two entries:
\itemize{
\item \strong{$parameter} The estimated parameters.
\item \strong{$mu} The estimated mean.
}
}
\description{
This function performs the estimation of the mean of a set of curves using
local linear smoothers where the bandwidth is estimated using the
methodology from Golovkine et al. (2021).
}
\references{
Golovkine S., Klutchnikoff N., Patilea V. (2021) - Adaptive
estimation of irregular mean and covariance functions.
}
