# Selection regimes

A simulation runs through a list of regimes:

## Usage

``` r
drift(generations)

stabilizing(generations, omega, theta = NULL)

directional(generations, delta_theta, omega, theta = NULL)

fluctuating(
  generations,
  every,
  step,
  omega,
  structure = "correlated",
  bounds = c(-1, 1),
  theta = NULL
)
```

## Arguments

- generations:

  number of generations.

- omega:

  fitness surface covariance matrix (traits x traits), or a single
  number `w`, meaning `w` times the identity.

- theta:

  optimum at the start of the regime, or `NULL`.

- delta_theta:

  change in the optimum per generation (vector, one value per trait).

- every:

  number of generations between jumps of the optimum.

- step:

  size of a jump (one value, or one per trait).

- structure:

  `"correlated"`, `"independent"` or a list of trait blocks.

- bounds:

  lower and upper bound of the optimum (`NULL` for none).

## Value

An `evomod_regime` object.

## Details

- `drift()`: no selection.

- `stabilizing()`: Gaussian fitness surface with optimum `theta` and
  covariance `omega`,
  `w(z) = exp(-(z - theta)' omega^-1 (z - theta) / 2)`.

- `directional()`: as stabilizing, but the optimum moves by
  `delta_theta` every generation (before selection).

- `fluctuating()`: as stabilizing, but after every `every` generations
  the optimum jumps by `step` up or down (do O & Whitlock 2023). With
  `structure = "correlated"` all traits jump in the same direction; with
  `"independent"` each trait draws its own direction; a list of trait
  blocks (e.g. `list(1:2, 3:4)`) makes traits jump together within
  blocks and independently between them (traits not in any block jump
  alone). A jump that would take a trait of the block outside `bounds`
  is not taken, so at a bound the optimum stays or moves inward with
  probability 1/2 each; `bounds = NULL` leaves the optimum unbounded.

When `theta` is `NULL` the regime starts from the optimum at the end of
the previous regime (the origin for the first one).

## Examples

``` r
regimes <- list(drift(100), stabilizing(1000, omega = 20),
                directional(1000, delta_theta = c(0.005, 0), omega = 20))
fluctuating(5000, every = 10, step = 0.2, omega = 1, structure = "correlated")
#> <evomod regime> fluctuating for 5000 generations
fluctuating(5000, every = 10, step = 0.2, omega = 1, structure = list(1:2, 3:4))
#> <evomod regime> fluctuating for 5000 generations
```
