# Recorded population statistics

Statistics recorded during
[`evolve()`](https://diogro.github.io/evomod2/reference/evolve.md),
every `record_every` generations.

## Usage

``` r
trait_means(pop)

mean_fitness(pop)

variances(pop)

g_matrices(pop, component = c("total", "genic", "ld"))

p_matrices(pop)
```

## Arguments

- pop:

  an `evomod_population`.

- component:

  `"total"` (G), `"genic"` or `"ld"`.

## Value

A data frame or a named list of matrices.

## Details

- `trait_means()`: mean genotypic value (`mean_x`), mean phenotype
  (`mean_z`) and optimum (`theta`) of each trait, in long format.

- `mean_fitness()`: mean fitness of the population on the fitness
  surface of the current regime (`NA` under drift).

- `variances()`: additive genetic (`G`) and phenotypic (`P`) variances
  and heritability (`h2`) of each trait, in long format.

- `g_matrices()`, `p_matrices()`: lists of covariance matrices of
  genotypic values (G) and phenotypes (P), named by generation. For
  architectures where genotypic values are sums of allele effects
  ([`additive()`](https://diogro.github.io/evomod2/reference/architectures.md),
  [`binary_b()`](https://diogro.github.io/evomod2/reference/architectures.md))
  G can be split into its genic part, the G expected without linkage
  disequilibrium (twice the sum over loci of the covariance of allele
  effects), and the rest, due to linkage disequilibrium between loci and
  departures from Hardy-Weinberg (`component = "ld"`); see also
  [`g_decomposition()`](https://diogro.github.io/evomod2/reference/g_decomposition.md).
