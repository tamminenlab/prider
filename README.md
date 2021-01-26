# Prider - primer designer for complex DNA sequence sets

Prider detects groups of conserved sequence areas across complex, internally diverse DNA sequence sets. Prider leverages the fast, scalable sequence similarity search algorithm implemented in [Blaster](https://github.com/manutamminen/blaster).

## Installation

```R
# install.packages("devtools")
devtools::install_github("manutamminen/prider")

```

## Examples

```R
primer_designs <- prider("test.fasta", primer_length = 20)
```
