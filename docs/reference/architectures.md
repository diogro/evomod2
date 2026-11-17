# Genetic architectures

An architecture maps a genome to genotypic values. Genomes are matrices
of allelic effect vectors (traits x alleles); alleles `2l - 1` and `2l`
belong to locus `l`, and the effect of a locus is the sum of its two
alleles (`Y_l`). Effects are measured from the reference genotype, where
all effects are zero.

## Usage

``` r
additive(loci, traits)

binary_b(loci, traits, init_b = 0.25, init_y = c("normal", "mutations"))

contribution_matrix(
  traits,
  genes = traits,
  init_weights = rep(1, genes),
  init_contributions = diag(1, traits, genes)
)

multilinear(loci, traits, epsilon)
```

## Arguments

- loci:

  number of diploid loci.

- traits:

  number of traits.

- init_b:

  probability that an entry of `B` is 1 in a new population (0.25 in the
  C code, which starts from random mutations of an empty `B`).

- init_y:

  how allelic effects `y` of a new population are drawn, with the
  model's `init_sd`: `"normal"`, `Normal(0, init_sd)` for every allele;
  `"mutations"`, as in the C code, `2 x loci` mutations of
  `Normal(0, init_sd)` added to alleles chosen at random in each genome.

- genes:

  number of weight genes.

- init_weights:

  starting values of the weight genes (one per gene).

- init_contributions:

  starting contribution matrix (`traits x genes`); by default ones on
  the diagonal (gene `j` affects trait `j` only, the identity when
  `genes == traits`).

- epsilon:

  epistasis coefficients: a `traits x traits x traits` array shared by
  all pairs of loci (symmetric in its last two indices; see
  [`random_epsilon()`](https://diogro.github.io/evomod2/reference/random_epsilon.md)),
  or a `traits x traits x traits x loci x loci` array with a different
  tensor for every pair of loci, as in Jones, Burger & Arnold (2014)
  (see
  [`pairwise_epsilon()`](https://diogro.github.io/evomod2/reference/pairwise_epsilon.md)).
  For the pairwise form, `epsilon[k, a, b, i, j]` multiplies
  `Y_i[a] Y_j[b]` in trait `k`, and only entries with `i < j` are used.

## Value

An `evomod_architecture` object.

## Details

- `additive()`: `x = sum_l Y_l`.

- `binary_b()`: the binary ontogenetic matrix model of Melo & Marroig
  (2015). Each allele has a scalar effect `y` and a column of a binary
  matrix `B` saying which traits it affects; `x = B y`. Use with
  [`binary_b_mutation()`](https://diogro.github.io/evomod2/reference/mutation.md).
  Genomes have `traits + 1` rows: `y`, then `B`.

- `contribution_matrix()`: the model of do O & Whitlock (2023). `genes`
  weight genes `a_j` and a `traits x genes` matrix of contribution genes
  `c_ij`, each its own locus with a scalar value (the two alleles of a
  locus are summed); `x_i = sum_j c_ij a_j`. Pleiotropy (`C`) evolves
  continuously and can be negative, and the map is bilinear (epistasis
  between weight and contribution genes). Loci `1..genes` are the weight
  genes; locus `genes + (j - 1) * traits + i` is `c_ij`. Use with
  [`contribution_mutation()`](https://diogro.github.io/evomod2/reference/mutation.md).
  New populations start with every individual at `init_weights` and
  `init_contributions` (each allele half of the locus value), plus
  `Normal(0, init_sd)` noise on every allele.

- `multilinear()` with a pairwise `epsilon`
  ([`pairwise_epsilon()`](https://diogro.github.io/evomod2/reference/pairwise_epsilon.md)):
  a different tensor for every pair of loci (Jones, Burger & Arnold
  2014), `x = S + sum_{i<j} eps_ij(Y_i, Y_j)`.

- `multilinear()`: pairwise multilinear epistasis (Hansen & Wagner
  2001), `x = S + sum_{i<j} eps(Y_i, Y_j)`, with `S = sum_l Y_l` and
  `eps(a, b)[k] = sum_{l,m} epsilon[k, l, m] a[l] b[m]`. The same tensor
  is used for every pair of loci and must be symmetric in its last two
  indices.

## Examples

``` r
additive(loci = 50, traits = 2)
#> <evomod architecture> additive - 50 loci, 2 traits
binary_b(loci = 50, traits = 4)
#> <evomod architecture> binary_b - 50 loci, 4 traits
contribution_matrix(traits = 4)
#> <evomod architecture> contribution_matrix - 20 loci, 4 traits
multilinear(loci = 50, traits = 2, epsilon = random_epsilon(2, sd = 0.05))
#> <evomod architecture> multilinear - 50 loci, 2 traits
```
