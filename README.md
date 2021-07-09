# Prider

Prider detects groups of conserved sequence areas across complex, internally diverse DNA sequence sets.

## Installation

```R
# install.packages("devtools")
devtools::install_github("manutamminen/prider")

```

## Examples

```R
test_fasta <- system.file("extdata", "test.fasta", package = "prider")

primer_designs <- prider(test_fasta)

primers(primer_designs)

primers(primer_designs)[1]

sequences(primer_designs)

sequences(primer_designs)[1]
```
