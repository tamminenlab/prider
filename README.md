# Prider

Prider detects groups of conserved sequence areas across complex, internally 
diverse DNA sequence sets.

## Installation

```R
# install.packages("devtools")
devtools::install_github("manutamminen/prider")

```

## Examples

```R
test_fasta <- system.file("extdata", "test.fasta", package = "prider")

# Runs Prider with default values:
primer_designs <- prider(test_fasta)

# Returns all the primers:
primers(primer_designs)

# Returns the primers of a specific primer group:
primers(primer_designs)[1]

# Returns all the sequences:
sequences(primer_designs)

# Returns the sequence of a specific Id:
sequences(primer_designs)[1]

# Plots the primers groups and the target sequences as a heatmap:
plot(primer_designs)

```
