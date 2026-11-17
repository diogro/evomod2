# Mutation models

- `mutation()`: for
  [`additive()`](https://diogro.github.io/evomod2/reference/architectures.md)
  and
  [`multilinear()`](https://diogro.github.io/evomod2/reference/architectures.md)
  architectures; `hoc_mutation()` is the house of cards version
  (`persistence = 0`). Each allele mutates with probability `rate` per
  generation; a mutation adds a draw from a multivariate normal with
  mean zero and covariance `cov` to the allele's effect vector.

- `contribution_mutation()`: for
  [`contribution_matrix()`](https://diogro.github.io/evomod2/reference/architectures.md).
  Each allele of a weight gene mutates with probability `rate_weight`,
  each allele of a contribution gene with `rate_contribution`; a
  mutation adds `Normal(0, sd)` to the allele's value.

- `binary_b_mutation()`: for
  [`binary_b()`](https://diogro.github.io/evomod2/reference/architectures.md).
  Each allele's effect `y` mutates with probability `rate`, adding
  `Normal(0, sd)`; each entry of `B` mutates with probability `rate_b`.
  With `b_mutation = "flip"` (the default, and the current C code) a
  mutated entry always changes; with `"random"` it gets a random value,
  0 or 1, so it changes half of the time, as in the C code used for the
  simulations of Melo & Marroig (2015).

## Usage

``` r
mutation(rate, cov, persistence = 1)

hoc_mutation(rate, cov)

contribution_mutation(rate_weight, rate_contribution, sd, persistence = 1)

binary_b_mutation(
  rate,
  sd,
  rate_b,
  b_mutation = c("flip", "random"),
  persistence = 1
)
```

## Arguments

- rate:

  mutation probability per allele per generation.

- cov:

  mutational covariance matrix (traits x traits), or a single number
  \\s\\, meaning \\s^2 I\\.

- persistence:

  how much of the old effect a mutation keeps: the new effect is \\\rho
  e + \delta\\, with \\\rho\\ = `persistence`, \\e\\ the old effect and
  \\\delta\\ a fresh draw. 1 (the default) is a random walk, where
  mutational effects accumulate; 0 is the house of cards model (Kingman
  1978; Turelli 1984), where a mutation replaces the effect with a fresh
  draw unrelated to it; values in between give a regression towards the
  reference genotype (Zeng & Cockerham 1993). For
  [`binary_b()`](https://diogro.github.io/evomod2/reference/architectures.md)
  it applies to `y`.

- rate_weight, rate_contribution:

  mutation probability per allele per generation of weight genes and of
  contribution genes.

- sd:

  standard deviation of mutations in `y`.

- rate_b:

  mutation probability per entry of `B` per generation.

- b_mutation:

  `"random"` or `"flip"`, see above.

## Value

An `evomod_mutation` object.

## Examples

``` r
mutation(rate = 0.005, cov = 0.05)
#> $type
#> [1] "gaussian"
#> 
#> $rate
#> [1] 0.005
#> 
#> $cov
#> [1] 0.05
#> 
#> $persistence
#> [1] 1
#> 
#> attr(,"class")
#> [1] "evomod_mutation"
hoc_mutation(rate = 0.005, cov = 0.5)   # house of cards
#> $type
#> [1] "gaussian"
#> 
#> $rate
#> [1] 0.005
#> 
#> $cov
#> [1] 0.5
#> 
#> $persistence
#> [1] 0
#> 
#> attr(,"class")
#> [1] "evomod_mutation"
binary_b_mutation(rate = 0.0005, sd = 0.2, rate_b = 0.0001)
#> $type
#> [1] "binary_b"
#> 
#> $rate
#> [1] 5e-04
#> 
#> $sd
#> [1] 0.2
#> 
#> $rate_b
#> [1] 1e-04
#> 
#> $b_mutation
#> [1] "flip"
#> 
#> $random_state
#> [1] FALSE
#> 
#> $persistence
#> [1] 1
#> 
#> attr(,"class")
#> [1] "evomod_mutation"
mutation(rate = 0.005, cov = matrix(c(1, 0.5, 0.5, 1), 2) * 0.01)
#> $type
#> [1] "gaussian"
#> 
#> $rate
#> [1] 0.005
#> 
#> $cov
#>       [,1]  [,2]
#> [1,] 0.010 0.005
#> [2,] 0.005 0.010
#> 
#> $persistence
#> [1] 1
#> 
#> attr(,"class")
#> [1] "evomod_mutation"
```
