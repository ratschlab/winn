% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/winn.R
\name{autocorrelation_correct}
\alias{autocorrelation_correct}
\title{Correct for Drift in Data Using Autocorrelation Correction}
\usage{
autocorrelation_correct(
  data,
  run_order = NULL,
  batch,
  lag = 20,
  test = "Ljung-Box",
  detrend = "mean",
  fdr_threshold = 0.05
)
}
\arguments{
\item{data}{A numeric matrix with rows representing metabolites and columns representing samples.}

\item{run_order}{An optional numeric vector representing the run order of the samples.}

\item{batch}{A numeric vector indicating the batch (or segment) assignment for each sample.}

\item{lag}{An integer specifying the lag to be used in the autocorrelation test.}

\item{test}{A character string specifying the autocorrelation test to use ("Ljung-Box" or "DW").}

\item{detrend}{A character string indicating the method for detrending ("mean" or "spline").}

\item{fdr_threshold}{A numeric value specifying the FDR threshold for significance.}
}
\value{
A numeric matrix with drift corrected.
}
\description{
This function corrects for drift effects in metabolomics data by detrending based on run order within each batch segment.
}
\examples{
your_data_matrix <- matrix(rnorm(200, mean = 100, sd = 15), nrow = 20)
batch <- rep(1:4, length.out = ncol(your_data_matrix))
run_order <- seq_len(ncol(your_data_matrix))
drift_corrected <- autocorrelation_correct(your_data_matrix, run_order, batch)
}
