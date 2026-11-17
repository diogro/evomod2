# Package index

## Models

Genetic architectures, mutation, recombination and life cycles.

- [`evomod_model()`](https://diogro.github.io/evomod2/reference/evomod_model.md)
  : Simulation model
- [`additive()`](https://diogro.github.io/evomod2/reference/architectures.md)
  [`binary_b()`](https://diogro.github.io/evomod2/reference/architectures.md)
  [`contribution_matrix()`](https://diogro.github.io/evomod2/reference/architectures.md)
  [`multilinear()`](https://diogro.github.io/evomod2/reference/architectures.md)
  : Genetic architectures
- [`random_epsilon()`](https://diogro.github.io/evomod2/reference/random_epsilon.md)
  : Random epistasis tensor
- [`pairwise_epsilon()`](https://diogro.github.io/evomod2/reference/pairwise_epsilon.md)
  : Random pairwise epistasis coefficients
- [`mutation()`](https://diogro.github.io/evomod2/reference/mutation.md)
  [`hoc_mutation()`](https://diogro.github.io/evomod2/reference/mutation.md)
  [`contribution_mutation()`](https://diogro.github.io/evomod2/reference/mutation.md)
  [`binary_b_mutation()`](https://diogro.github.io/evomod2/reference/mutation.md)
  : Mutation models
- [`free_recombination()`](https://diogro.github.io/evomod2/reference/recombination.md)
  [`genetic_map()`](https://diogro.github.io/evomod2/reference/recombination.md)
  [`linkage()`](https://diogro.github.io/evomod2/reference/recombination.md)
  : Recombination
- [`wright_fisher()`](https://diogro.github.io/evomod2/reference/life_cycles.md)
  [`moran()`](https://diogro.github.io/evomod2/reference/life_cycles.md)
  [`viability_selection()`](https://diogro.github.io/evomod2/reference/life_cycles.md)
  : Life cycles

## Selection regimes and simulation

- [`drift()`](https://diogro.github.io/evomod2/reference/regimes.md)
  [`stabilizing()`](https://diogro.github.io/evomod2/reference/regimes.md)
  [`directional()`](https://diogro.github.io/evomod2/reference/regimes.md)
  [`fluctuating()`](https://diogro.github.io/evomod2/reference/regimes.md)
  : Selection regimes
- [`evolve()`](https://diogro.github.io/evomod2/reference/evolve.md) :
  Evolve a population
- [`plot(`*`<evomod_population>`*`)`](https://diogro.github.io/evomod2/reference/plot.evomod_population.md)
  : Plot recorded statistics

## Recorded statistics and genomes

- [`trait_means()`](https://diogro.github.io/evomod2/reference/statistics.md)
  [`mean_fitness()`](https://diogro.github.io/evomod2/reference/statistics.md)
  [`variances()`](https://diogro.github.io/evomod2/reference/statistics.md)
  [`g_matrices()`](https://diogro.github.io/evomod2/reference/statistics.md)
  [`p_matrices()`](https://diogro.github.io/evomod2/reference/statistics.md)
  : Recorded population statistics
- [`genomes()`](https://diogro.github.io/evomod2/reference/genomes.md)
  [`genotypic_values()`](https://diogro.github.io/evomod2/reference/genomes.md)
  [`local_effects()`](https://diogro.github.io/evomod2/reference/genomes.md)
  : Genomes and genotypic values

## Analysis

G matrix decomposition, mutational covariance and matrix statistics.

- [`g_decomposition()`](https://diogro.github.io/evomod2/reference/g_decomposition.md)
  : Pleiotropic and linkage disequilibrium parts of G
- [`shuffle_genomes()`](https://diogro.github.io/evomod2/reference/shuffle_genomes.md)
  : Shuffle alleles among individuals
- [`mutational_covariance()`](https://diogro.github.io/evomod2/reference/mutational_covariance.md)
  : Mutational covariance
- [`additive_covariance()`](https://diogro.github.io/evomod2/reference/additive_covariance.md)
  : Additive genetic covariance
- [`avg_corr()`](https://diogro.github.io/evomod2/reference/avg_corr.md)
  : Average correlations within and between modules
- [`autonomy()`](https://diogro.github.io/evomod2/reference/matrix_statistics.md)
  [`mean_autonomy()`](https://diogro.github.io/evomod2/reference/matrix_statistics.md)
  [`eigen_variance()`](https://diogro.github.io/evomod2/reference/matrix_statistics.md)
  : Matrix statistics: autonomy and eigenvalue distribution

## Package

- [`evomod2`](https://diogro.github.io/evomod2/reference/evomod2-package.md)
  [`evomod2-package`](https://diogro.github.io/evomod2/reference/evomod2-package.md)
  : evomod2: individual-based simulations of genetic architecture
  evolution
