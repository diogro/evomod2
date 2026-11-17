# Life cycles

- `wright_fisher()` (the default): discrete generations. Every
  generation, each individual mutates and the whole population is
  replaced by `n` offspring of parents drawn with probability
  proportional to fitness.

- `moran()`: overlapping generations. A generation is `n` birth-death
  events; at each event a random individual dies and is replaced by the
  offspring of two parents drawn with probability proportional to
  fitness (they can be the same individual, as in `wright_fisher()`).
  The offspring mutates with the model's per-generation rates, so the
  mutational input per generation is the same as with `wright_fisher()`.
  Drift per generation is about twice as strong as with
  `wright_fisher()` at the same `n`.

## Usage

``` r
wright_fisher()

moran()

viability_selection(mates = 2, offspring_per_mate = 2)
```

## Arguments

- mates:

  number of males each female mates with.

- offspring_per_mate:

  number of offspring per female and male.

## Value

An `evomod_life_cycle` object, for
[`evomod_model()`](https://diogro.github.io/evomod2/reference/evomod_model.md).

## Details

- `viability_selection()`: the life cycle of Jones, Arnold & Burger
  (2003, 2004, 2014). Half the population is female; each female mates
  with `mates` males drawn at random (males can mate many times) and has
  `offspring_per_mate` offspring with each. Zygotes mutate, survive with
  probability equal to their fitness (absolute, `W(z) <= 1`), and `n`
  adults are drawn at random from the survivors (if fewer than `n`
  survive, survivors are resampled to keep `n`). `n` must be even. With
  the defaults (2 mates, 2 offspring each), the effective population
  size is about `n`, as with `wright_fisher()`.

Statistics are recorded, and directional or fluctuating optima move,
between generations.

## Examples

``` r
moran()
#> $type
#> [1] "moran"
#> 
#> attr(,"class")
#> [1] "evomod_life_cycle"
```
