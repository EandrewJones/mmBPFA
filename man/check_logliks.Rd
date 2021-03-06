% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_results.R
\name{check_logliks}
\alias{check_logliks}
\title{Function to examine posterior log-Likelihoods}
\usage{
check_logliks(
  log_liks,
  probs = c(0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95),
  plot = TRUE
)
}
\arguments{
\item{log_liks}{A numeric vector of posterior draws from mmBPFA sampler. Must not be processed.}

\item{probs}{A numeric vector of quantiles to display for log-likelihoods distribution. Default is c(0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95).}

\item{plot}{A boolean indicating whether to plot the log-likelihood distribution. Default is TRUE.}
}
\value{
If plot is TRUE, returns none, else returns mean log-likelihood and quantiles.
}
\description{
\code{check_logliks} ingests the raw vector of posterior log-likelihoods and returns quantiles of interests, the mean, and
a plot of the distribution if \code{plot = TRUE}. A well-behaving posterior distribution should be unimodal and peak at the mean.
}
