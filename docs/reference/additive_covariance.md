# Additive genetic covariance

Least-squares additive genetic (co)variance of the current population:
the covariance matrix of the best linear prediction of genotypic values
from the state of every locus (the sum of its two alleles: their effect
vectors, or for
[`binary_b()`](https://diogro.github.io/evomod2/reference/architectures.md)
their effects `B y`, or for
[`contribution_matrix()`](https://diogro.github.io/evomod2/reference/architectures.md)
the gene values). For the architectures in evomod2 the genotypic value
is linear in the state of each locus given the others, so this is the
additive variance in the usual sense (variance of breeding values) under
linkage equilibrium. `G - additive` is the non-additive (epistatic)
covariance. Needs more individuals than loci times traits.

## Usage

``` r
additive_covariance(pop)
```

## Arguments

- pop:

  an `evomod_population`.

## Value

A list with `G`, `additive` and `nonadditive` (`G - additive`), all with
the `n - 1` denominator of
[`stats::cov()`](https://rdrr.io/r/stats/cor.html).

## Examples

``` r
eps <- pairwise_epsilon(loci = 10, traits = 2, sd = 1)
m <- evomod_model(multilinear(10, 2, eps), mutation(0.01, 0.2), env_sd = 1, init_sd = 0.1)
pop <- evolve(m, stabilizing(200, omega = 20), n = 300, seed = 1)
additive_covariance(pop)
#> $G
#>           trait1    trait2
#> trait1 1.8910878 0.1247169
#> trait2 0.1247169 1.6747977
#> 
#> $additive
#>           trait1    trait2
#> trait1 0.9021293 0.1499882
#> trait2 0.1499882 0.8469063
#> 
#> $nonadditive
#>             trait1      trait2
#> trait1  0.98895850 -0.02527132
#> trait2 -0.02527132  0.82789140
#> 
```
