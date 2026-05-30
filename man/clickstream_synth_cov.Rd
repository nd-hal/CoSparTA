\name{clickstream_synth_cov}
\alias{clickstream_synth_cov}
\docType{data}
\title{Synthetic clickstream covariate matrix}
\description{
A data frame of demographic covariates for the 4000-user synthetic
clickstream dataset bundled with CoSparTA. Each row corresponds to one
user. Variable distributions mirror the comScore panel structure used to
generate the synthetic tensor.
}
\format{A data frame with 4000 rows and 6 columns:
\describe{
  \item{gender_binary}{Integer. Binarized gender indicator (0/1).}
  \item{age_centered}{Numeric. Age centered and standardized.}
  \item{race_binary}{Integer. Binarized race indicator (0/1).}
  \item{hh_income_num}{Integer. Household income category (ordinal).}
  \item{hh_edu_num}{Integer. Household education level (ordinal).}
  \item{children_binary}{Integer. Presence of children in household (0/1/2).}
}}
\source{Simulated data generated to mirror comScore panel demographics.}
\keyword{datasets}
