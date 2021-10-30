  <!-- badges: start -->
  [![R-CMD-check](https://github.com/manutamminen/prider/workflows/R-CMD-check/badge.svg)](https://github.com/manutamminen/prider/actions)
  [![CRAN status](https://www.r-pkg.org/badges/version/prider)](https://CRAN.R-project.org/package=prider)
  <!-- badges: end -->

# Prider

Prider permits multiplexed oligonucleotide primer and probe design for 
complex DNA sequence sets by implementing an algorithm for linearly 
scaling approximation of set coverage. A detailed description available at
[Smolander and Tamminen, 2021](https://www.biorxiv.org/content/10.1101/2021.09.06.459073v1).

## Installation

```R
# install.packages("devtools")
devtools::install_github("manutamminen/prider")

```

## Examples

```R
test_fasta <- system.file("extdata", "test.fasta", package = "prider")

# Runs Prider with default values
primer_designs <- prider(test_fasta)

# Returns all the primers
primers(primer_designs)

# Returns the primers of a specific primer group
primers(primer_designs)[1]

# Returns all the sequences
sequences(primer_designs)

# Returns the sequence of a specific Id
sequences(primer_designs)[1]

# Plots the primers groups and the target sequences as a heatmap
plot(primer_designs)

```
