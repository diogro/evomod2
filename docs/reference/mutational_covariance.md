# Mutational covariance

Effects of new mutations on genotypic values in the current population.
Genomes are drawn at random (with replacement), and each receives
exactly one mutation event of the model's mutation process, at a site
chosen in proportion to the mutation rates; the effect is the change in
its genotypic value. With epistasis
([`multilinear()`](https://diogro.github.io/evomod2/reference/architectures.md))
or
[`binary_b()`](https://diogro.github.io/evomod2/reference/architectures.md)
the effects depend on the genetic background, so the mutational
covariance evolves.

## Usage

``` r
mutational_covariance(pop, n = 10000, effects = FALSE)
```

## Arguments

- pop:

  an `evomod_population`.

- n:

  number of mutations sampled.

- effects:

  if `TRUE`, also return the sampled effects.

## Value

A list with

- `mean`: mean effect of a mutation on each trait (non-zero with
  directional epistasis);

- `per_mutation`: covariance matrix of the effects of one mutation;

- `expected_mutations`: expected number of mutation events per genome
  per generation;

- `M`: mutational (co)variance input per generation, \\\lambda\\
  E\[\Delta x \Delta x^\top\]\\, with \\\lambda\\ =
  `expected_mutations`;

- `effects`: the sampled effects (traits x n), if requested.

## Examples

``` r
model <- evomod_model(additive(loci = 20, traits = 2),
                      mutation(rate = 0.01, cov = 0.1), init_sd = 0.1)
pop <- evolve(model, stabilizing(20, omega = 10), n = 100, seed = 1)
mutational_covariance(pop, n = 5000)
#> $mean
#>        trait1        trait2 
#>  0.0008153449 -0.0021901519 
#> 
#> $per_mutation
#>              trait1       trait2
#> trait1 9.639040e-03 5.277904e-05
#> trait2 5.277904e-05 1.007471e-02
#> 
#> $expected_mutations
#> [1] 0.4
#> 
#> $M
#>              trait1       trait2
#> trait1 0.0038551108 0.0000203931
#> trait2 0.0000203931 0.0040309950
#> 
```
