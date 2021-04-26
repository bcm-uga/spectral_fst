# spectraltheory

## Summary

A package to compute values and matrices of interest described in the paper "A Spectral Theory for Wright’s Inbreeding Coefficients and Related Quantities" François et Gain 2020. It allows you to compute Fst of surrogate genotype. For instance, you can obtain the Fst after removing adaptive genetic variation associated with environmental variables.

## Package overview

The main function of the package is

* **compute_partition** Computing partition matrices and related quantities from spectral theory


## Installation

Install the latest version from github (requires [devtools](https://github.com/hadley/devtools)):
```R
# install.packages("devtools")
devtools::install_github("bcm-uga/spectralfst")
```

## References

- Olivier François, Clément Gain. (2021). A Spectral Theory for Wright’s Inbreeding Coefficients and Related Quantities. In review.
