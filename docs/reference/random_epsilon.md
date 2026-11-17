# Random epistasis tensor

Entries are drawn from `Normal(mean, sd)` and symmetrised in the last
two indices. A positive `mean` gives positive directional epistasis on
average.

## Usage

``` r
random_epsilon(traits, sd = 1, mean = 0)
```

## Arguments

- traits:

  number of traits.

- sd, mean:

  standard deviation and mean of the entries.

## Value

A `traits x traits x traits` array.

## Examples

``` r
random_epsilon(2, sd = 0.1)
#> , , 1
#> 
#>             [,1]       [,2]
#> [1,]  0.05416082 0.10730829
#> [2,] -0.09958599 0.09058614
#> 
#> , , 2
#> 
#>            [,1]       [,2]
#> [1,] 0.10730829 -0.0734806
#> [2,] 0.09058614  0.1322102
#> 
```
