% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/winn.R
\name{normalize_by_dilution_factor}
\alias{normalize_by_dilution_factor}
\title{Normalize Data Matrix by Dilution Factor}
\usage{
normalize_by_dilution_factor(
  data,
  processing = "shrink",
  control_samples = NULL
)
}
\arguments{
\item{data}{A numeric matrix or data frame where rows represent metabolites and columns represent samples.}

\item{processing}{A character string specifying the processing method to use. Options are "shrink" (default) or "normalize".}

\item{control_samples}{An optional numeric vector specifying which columns correspond to control samples.
If provided, the reference spectrum is calculated using only these samples.}
}
\value{
A numeric matrix of normalized intensities.
}
\description{
This function normalizes a data matrix by dilution factors or alternatively shrinks sample measurements
if any samples have dilution factors that are more than one standard deviation away from the mean dilution factor.
}
\examples{
# Example usage:
your_data_matrix <- matrix(rnorm(200, mean = 100, sd = 15), nrow = 20)
normalized_data <- normalize_by_dilution_factor(your_data_matrix, control_samples = 1:4)
}
