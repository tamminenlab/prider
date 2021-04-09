# Prider

Prider detects groups of conserved sequence areas across complex, internally diverse DNA sequence sets.

## Installation

```R
# install.packages("devtools")
devtools::install_github("manutamminen/prider")

```

## Examples

```R
primer_designs <- prider("test.fasta", primer_length = 20)
```
