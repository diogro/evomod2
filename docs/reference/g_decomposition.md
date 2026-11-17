# Pleiotropic and linkage disequilibrium parts of G

Splits the G matrix of the current population into the G expected
without associations between loci (`genic`) and the remainder (`ld`),
due to linkage disequilibrium between loci and departures from
Hardy-Weinberg. All matrices use the `n - 1` denominator of
[`stats::cov()`](https://rdrr.io/r/stats/cor.html), so `G = genic + ld`.

## Usage

``` r
g_decomposition(pop, method = c("genic", "shuffle"), shuffles = 10)
```

## Arguments

- pop:

  an `evomod_population`.

- method:

  `"genic"` or `"shuffle"`, see above.

- shuffles:

  number of shuffled populations averaged (`method = "shuffle"`).

## Value

A list with matrices `G`, `genic` and `ld`.

## Details

- `method = "genic"` (exact; architectures where genotypic values are
  sums of allele effects,
  [`additive()`](https://diogro.github.io/evomod2/reference/architectures.md)
  and
  [`binary_b()`](https://diogro.github.io/evomod2/reference/architectures.md)):
  twice the sum over loci of the covariance of allele effects, pooling
  both alleles of every individual.

- `method = "shuffle"` (any architecture, including
  [`multilinear()`](https://diogro.github.io/evomod2/reference/architectures.md)):
  at each locus the `2n` allele copies are pooled and randomly
  reassigned to individuals, which removes associations between loci and
  within loci but keeps allele frequencies and every allele's full state
  (its effect vector, or its `y` and column of `B`). `genic` is the mean
  G of `shuffles` shuffled populations (do O & Whitlock 2023 call it the
  pleiotropic part of G). For additive architectures its expectation is
  the exact genic part.

## References

do O, I., & Whitlock, M. C. (2023). The evolution of genetic covariance
and modularity as a result of multigenerational environmental
fluctuation. *Evolution Letters*, 7(6), 457-466.

## Examples

``` r
model <- evomod_model(additive(loci = 20, traits = 2),
                      mutation(rate = 0.01, cov = 0.1), init_sd = 0.1)
pop <- evolve(model, stabilizing(20, omega = 10), n = 100, seed = 1)
g_decomposition(pop)
#> $G
#>            trait1     trait2
#> trait1 0.47362851 0.02699102
#> trait2 0.02699102 0.40358214
#> 
#> $genic
#>             trait1      trait2
#> trait1  0.51782676 -0.02774668
#> trait2 -0.02774668  0.44735322
#> 
#> $ld
#>             trait1      trait2
#> trait1 -0.04419825  0.05473770
#> trait2  0.05473770 -0.04377108
#> 
g_decomposition(pop, method = "shuffle", shuffles = 20)
#> $G
#>            trait1     trait2
#> trait1 0.47362851 0.02699102
#> trait2 0.02699102 0.40358214
#> 
#> $genic
#>              trait1       trait2
#> trait1  0.530295375 -0.007434872
#> trait2 -0.007434872  0.452740317
#> 
#> $ld
#>             trait1      trait2
#> trait1 -0.05666687  0.03442589
#> trait2  0.03442589 -0.04915818
#> 
```
