# Recombination

How the alleles a parent passes on are put together:

## Usage

``` r
free_recombination()

genetic_map(chromosomes = 1, length = 1)

linkage(rates)
```

## Arguments

- chromosomes:

  number of chromosomes.

- length:

  length of each chromosome in Morgans (1 Morgan = 100 cM).

- rates:

  recombination fractions between loci `l` and `l + 1`.

## Value

An `evomod_recombination` object, resolved for the number of loci by
[`evomod_model()`](https://diogro.github.io/evomod2/reference/evomod_model.md).

## Details

- `free_recombination()` (the default): every locus is inherited
  independently.

- `genetic_map()`: loci are placed in order, evenly spaced, on
  `chromosomes` chromosomes of equal length (`length`, in Morgans; loci
  `1..k` on the first chromosome, and so on). Recombination fractions
  between neighbouring loci follow Haldane's map function, \\r = (1 -
  e^{-2d}) / 2\\ for a distance \\d\\; loci on different chromosomes
  recombine with probability 1/2.

- `linkage()`: explicit recombination fractions between neighbouring
  loci, a vector of length `loci - 1` (0.5 for unlinked loci).

Columns `2l - 1` and `2l` of a genome are the paternal and maternal
haplotypes at locus `l`; a gamete starts on a random haplotype and
switches at each crossover.

## Examples

``` r
genetic_map(chromosomes = 2, length = 0.5)
#> <evomod recombination> 2 chromosome(s) of 0.5 Morgans 
linkage(c(0.01, 0.5, 0.01))
#> <evomod recombination> recombination fractions for 4 loci 
```
