# Random pairwise epistasis coefficients

A different epistasis tensor for every pair of loci, with entries drawn
from `Normal(mean, sd)`, as in Jones, Burger & Arnold (2014), where the
variance of the coefficients is `sd^2`. With `within_trait = TRUE` only
the coefficients for interactions of a trait with itself
(`epsilon[k, k, k, i, j]`) are non-zero; `traits_with_epistasis`
restricts epistasis to some traits.

## Usage

``` r
pairwise_epsilon(
  loci,
  traits,
  sd = 1,
  mean = 0,
  within_trait = FALSE,
  traits_with_epistasis = seq_len(traits)
)
```

## Arguments

- loci:

  number of loci.

- traits:

  number of traits.

- sd, mean:

  standard deviation and mean of the coefficients.

- within_trait:

  only interactions of each trait with itself.

- traits_with_epistasis:

  traits whose values have epistatic terms.

## Value

A `traits x traits x traits x loci x loci` array (zero for `i >= j`),
for
[`multilinear()`](https://diogro.github.io/evomod2/reference/architectures.md).

## Examples

``` r
eps <- pairwise_epsilon(loci = 10, traits = 2, sd = 1)
multilinear(loci = 10, traits = 2, epsilon = eps)
#> <evomod architecture> multilinear_pairwise - 10 loci, 2 traits
```
