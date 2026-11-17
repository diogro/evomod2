# Evolve a population

Runs a population through a list of selection regimes. Starting from a
model creates a new population of `n` individuals; starting from a
population continues it (same individuals, generation count, optimum and
recorded statistics).

## Usage

``` r
evolve(x, regimes, ...)

# S3 method for class 'evomod_model'
evolve(x, regimes, n, record_every = 1, seed = NULL, initial = NULL, ...)

# S3 method for class 'evomod_population'
evolve(x, regimes, record_every = x$record_every, seed = NULL, ...)
```

## Arguments

- x:

  an `evomod_model` (new population) or an `evomod_population` (continue
  it).

- regimes:

  a regime or a list of regimes, see
  [regimes](https://diogro.github.io/evomod2/reference/regimes.md).

- ...:

  further arguments for the methods.

- n:

  population size.

- record_every:

  record population statistics every this many generations (and at the
  start).

- seed:

  optional seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the run.

- initial:

  optional starting genomes, an array `rows x alleles x n` (see
  [`genomes()`](https://diogro.github.io/evomod2/reference/genomes.md));
  by default they are drawn as described in
  [`evomod_model()`](https://diogro.github.io/evomod2/reference/evomod_model.md).

## Value

An `evomod_population`; see
[`trait_means()`](https://diogro.github.io/evomod2/reference/statistics.md),
[`g_matrices()`](https://diogro.github.io/evomod2/reference/statistics.md)
and
[`variances()`](https://diogro.github.io/evomod2/reference/statistics.md)
for the recorded statistics.

## Details

Each generation: the optimum moves (directional regimes), fitness is
computed from the current phenotypes, every individual mutates, `n`
pairs of parents are drawn with probability proportional to fitness
(with replacement), and each offspring gets one random allele per locus
from each parent (free recombination). Random numbers come from R's
generator, so [`set.seed()`](https://rdrr.io/r/base/Random.html) (or
`seed`) makes runs reproducible.

## Examples

``` r
model <- evomod_model(additive(loci = 20, traits = 2),
                      mutation(rate = 0.01, cov = 0.1), env_sd = 1, init_sd = 0.1)
pop <- evolve(model, n = 100, regimes = list(stabilizing(50, omega = 10)),
              record_every = 10, seed = 1)
pop <- evolve(pop, directional(50, delta_theta = c(0.02, 0), omega = 10))
trait_means(pop)
#>    generation  trait       mean_x      mean_z theta
#> 1           0 trait1 -0.105395381 -0.08994999   0.0
#> 2           0 trait2  0.052851082 -0.01485952   0.0
#> 3          10 trait1 -0.078026866 -0.06467254   0.0
#> 4          10 trait2  0.350370869  0.31476291   0.0
#> 5          20 trait1  0.001121080 -0.05145684   0.0
#> 6          20 trait2 -0.013594507  0.01421813   0.0
#> 7          30 trait1 -0.058305861 -0.06328598   0.0
#> 8          30 trait2 -0.323667988 -0.39672571   0.0
#> 9          40 trait1 -0.239821604 -0.11345083   0.0
#> 10         40 trait2 -0.099918544 -0.03504065   0.0
#> 11         50 trait1  0.286186814  0.43058418   0.0
#> 12         50 trait2 -0.214598556 -0.46836788   0.0
#> 13         60 trait1 -0.170151772 -0.26868587   0.2
#> 14         60 trait2 -0.072109369  0.17431115   0.0
#> 15         70 trait1  0.065662606  0.12602295   0.4
#> 16         70 trait2 -0.082971219 -0.10867104   0.0
#> 17         80 trait1  0.268709503  0.29918944   0.6
#> 18         80 trait2 -0.003612657 -0.05265973   0.0
#> 19         90 trait1  0.556115225  0.60307787   0.8
#> 20         90 trait2  0.250311010  0.23281721   0.0
#> 21        100 trait1  0.963869701  0.85990118   1.0
#> 22        100 trait2 -0.090043897 -0.23026369   0.0
```
