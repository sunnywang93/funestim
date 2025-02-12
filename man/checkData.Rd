% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{checkData}
\alias{checkData}
\title{Check the input data and return a list in the right format.}
\usage{
checkData(data)
}
\arguments{
\item{data}{An object from the package \code{funData}. It could be a
\code{funData::funData} or \code{funData::irregFunData} object.}
}
\value{
A list, where each element represents a curve. Each curve is defined
as a list with two entries:
\itemize{
\item \strong{$t} The sampling points
\item \strong{$x} The observed points
}
}
\description{
Check the input data and return a list in the right format.
}
\references{
C. Happ-Kurz (2020) Object-Oriented Software for Functional Data.
Journal of Statistical Software, 93(5): 1-38 .
}
