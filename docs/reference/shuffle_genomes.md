# Shuffle alleles among individuals

At each locus, pools the `2n` allele copies of the population and
assigns them to individuals at random. Allele frequencies and the state
of every allele are kept; associations between loci (linkage
disequilibrium) and between the two alleles of a locus (departures from
Hardy-Weinberg) are removed.

## Usage

``` r
shuffle_genomes(alleles)
```

## Arguments

- alleles:

  genomes, an array `rows x alleles x individuals` (see
  [`genomes()`](https://diogro.github.io/evomod2/reference/genomes.md)).

## Value

An array of the same shape.

## Examples

``` r
model <- evomod_model(additive(loci = 5, traits = 2), mutation(0.01, 0.1), init_sd = 0.1)
pop <- evolve(model, drift(5), n = 20, seed = 1)
shuffled <- shuffle_genomes(genomes(pop))
```
