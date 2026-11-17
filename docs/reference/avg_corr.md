# Average correlations within and between modules

Average correlations within and between modules

## Usage

``` r
avg_corr(m, modules = NULL)
```

## Arguments

- m:

  a correlation (or covariance) matrix.

- modules:

  list of trait indices of each module; by default the first and second
  half of the traits.

## Value

A list with the mean correlation within each module (`within`), the mean
correlation between traits in different modules (`between`), and their
ratio (`avg_ratio`).

## Examples

``` r
m <- matrix(c(1, .8, .1, .1, .8, 1, .1, .1, .1, .1, 1, .6, .1, .1, .6, 1), 4)
avg_corr(m)
#> $within
#> [1] 0.8 0.6
#> 
#> $between
#> [1] 0.1
#> 
#> $avg_ratio
#> [1] 7
#> 
```
