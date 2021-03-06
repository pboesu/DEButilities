% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_lb.R
\name{get_lb}
\alias{get_lb}
\title{get_lb}
\usage{
get_lb(p, eb = 1, lb0 = NA)
}
\arguments{
\item{p}{3-vector with parameters: g, k, v_H^b (see below)}

\item{eb}{optional scalar with scaled reserve density at birth (default eb = 1)}

\item{lb0}{optional scalar with initial estimate for scaled length at birth (default lb0 = NA, will use lb for k = 1)}
}
\value{
[lb, info] = lb: scalar with scaled length at birth; info: indicator equals 1 if successful, 0 otherwise
}
\description{
Obtains scaled length at birth, given the scaled reserve density at birth
}
\author{
Bas Kooijman
}
