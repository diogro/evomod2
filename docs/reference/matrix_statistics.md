# Matrix statistics: autonomy and eigenvalue distribution

- `autonomy()`: autonomy of a covariance matrix in the direction(s)
  `beta` (Hansen & Houle 2008),
  `a(beta) = 1 / ((beta' G^-1 beta) (beta' G beta))` for unit vectors
  `beta`: the fraction of the variation in direction `beta` that is
  independent of other directions (1 when `beta` is an eigenvector).

- `mean_autonomy()`: mean autonomy over `n` random directions.

- `eigen_variance()`: proportion of the total variance on each principal
  component.

## Usage

``` r
autonomy(m, beta)

mean_autonomy(m, n = 1000)

eigen_variance(m)
```

## Arguments

- m:

  a covariance matrix (e.g. a G or P matrix).

- beta:

  a direction (vector) or several directions (columns of a matrix);
  normalised to unit length.

- n:

  number of random directions.

## Value

`autonomy()`: one value per direction. `mean_autonomy()`: a number.
`eigen_variance()`: a vector of proportions in decreasing order.

## Examples

``` r
G <- matrix(c(1, 0.8, 0.1, 0.8, 1, 0.1, 0.1, 0.1, 1), 3)
autonomy(G, c(1, 1, 0))
#> [1] 0.9888889
mean_autonomy(G)
#> [1] 0.637667
eigen_variance(G)
#> [1] 0.60808802 0.32524531 0.06666667
```
