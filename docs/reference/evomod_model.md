# Simulation model

Combines a genetic architecture, a mutation model, recombination and
environmental noise.

## Usage

``` r
evomod_model(
  architecture,
  mutation,
  env_sd = 0,
  init_sd = 0,
  recombination = free_recombination(),
  life_cycle = wright_fisher()
)
```

## Arguments

- architecture:

  an architecture, see
  [architectures](https://diogro.github.io/evomod2/reference/architectures.md).

- mutation:

  a mutation model matching the architecture, see
  [`mutation()`](https://diogro.github.io/evomod2/reference/mutation.md).

- env_sd:

  standard deviation of environmental effects \\e\\ in \\z = x + e\\.

- init_sd:

  standard deviation of allelic effects in a new population (`y` for
  [`binary_b()`](https://diogro.github.io/evomod2/reference/architectures.md));
  0 (the default) starts every individual at the reference genotype.

- recombination:

  how loci are inherited, see
  [recombination](https://diogro.github.io/evomod2/reference/recombination.md);
  free recombination by default.

- life_cycle:

  [`wright_fisher()`](https://diogro.github.io/evomod2/reference/life_cycles.md)
  (discrete generations, the default),
  [`moran()`](https://diogro.github.io/evomod2/reference/life_cycles.md)
  (overlapping generations) or
  [`viability_selection()`](https://diogro.github.io/evomod2/reference/life_cycles.md)
  (separate sexes and viability selection, Jones et al.).

## Value

An `evomod_model` object.

## Examples

``` r
model <- evomod_model(additive(loci = 50, traits = 2),
                      mutation(rate = 0.005, cov = 0.05), env_sd = 1)
```
