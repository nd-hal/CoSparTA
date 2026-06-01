# Synthetic clickstream covariate matrix

A data frame of demographic covariates for the 4000-user synthetic
clickstream dataset bundled with CoSparTA. Each row corresponds to one
user. Variable distributions mirror the comScore panel structure used to
generate the synthetic tensor.

## Format

A data frame with 4000 rows and 6 columns:

- gender_binary:

  Integer. Binarized gender indicator (0/1).

- age_centered:

  Numeric. Age centered and standardized.

- race_binary:

  Integer. Binarized race indicator (0/1).

- hh_income_num:

  Integer. Household income category (ordinal).

- hh_edu_num:

  Integer. Household education level (ordinal).

- children_binary:

  Integer. Presence of children in household (0/1/2).

## Source

Simulated data generated to mirror comScore panel demographics.
