# Genomes and genotypic values

- `genomes()`: genomes of the current population, an array
  `rows x alleles x individuals` (alleles `2l - 1` and `2l` are locus
  `l`). Rows are the allelic effect vector (one per trait), or, for
  [`binary_b()`](https://diogro.github.io/evomod2/reference/architectures.md),
  `y` followed by the column of `B`.

- `genotypic_values()`: genotypic values of genomes under an
  architecture (traits x individuals).

- `local_effects()`: effect of each locus of one genome in the
  background of the rest of the genome (traits x loci): the change in
  genotypic value when the locus is set to the reference genotype. Equal
  to the locus effects for the additive architecture.

## Usage

``` r
genomes(pop)

genotypic_values(x, alleles = genomes(x))

local_effects(x, alleles = genome_of(x, individual), individual = 1)
```

## Arguments

- pop:

  an `evomod_population`.

- x:

  an architecture, model or population.

- alleles:

  a genome (rows x alleles) or an array of genomes (rows x alleles x
  individuals), see `genomes()`.

- individual:

  index of an individual of `pop`.
