# evomod2

<!-- badges: start -->
[![R-CMD-check](https://github.com/diogro/evomod2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/diogro/evomod2/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Individual-based simulations of the evolution of genetic architecture in
multivariate quantitative traits, with additive, multilinear epistatic or
binary ontogenetic matrix genotype-phenotype maps. The simulation engine is written in C++ (Rcpp +
RcppArmadillo); everything is driven from R.

evomod2 extends the model of Melo & Marroig (2015, *PNAS*
[10.1073/pnas.1322632112](https://doi.org/10.1073/pnas.1322632112)). The Julia
package [EvoMod.jl](https://github.com/diogro/evomod_julia) is the reference
implementation used for prototyping new models and for cross-checking this
package.

## Installation

```r
# install.packages("remotes")
remotes::install_github("diogro/evomod2")
```

A C++ compiler is needed (Xcode command line tools on macOS, Rtools on
Windows).

## Usage

A model combines a genetic architecture, a mutation model and environmental
noise. A population then evolves through a list of selection regimes.

```r
library(evomod2)

# pairwise multilinear epistasis between 50 loci affecting 2 traits
eps <- array(0, c(2, 2, 2))
eps[1, 1, 1] <- 0.05            # positive directional epistasis on trait 1
model <- evomod_model(
  architecture = multilinear(loci = 50, traits = 2, epsilon = eps),
  mutation = mutation(rate = 0.005, cov = 0.05),
  env_sd = 1
)

pop <- evolve(model, n = 1000, record_every = 20, seed = 1,
              regimes = list(stabilizing(2000, omega = 20),
                             directional(2000, delta_theta = c(0.005, 0), omega = 20)))

plot(pop)                 # mean genotypic values and optima
plot(pop, "variances")    # additive genetic variances
trait_means(pop)          # recorded statistics as data frames
variances(pop)
g_matrices(pop)           # G matrices by generation

pop <- evolve(pop, stabilizing(500, omega = 20))   # continue the same population
```

### Genic and linkage disequilibrium parts of G

For architectures where genotypic values are sums of allele effects
(`additive()`, `binary_b()`), G splits into a genic part, the G expected
without linkage disequilibrium (twice the sum over loci of the covariance of
allele effects), and the part due to linkage disequilibrium between loci:

```r
model <- evomod_model(additive(loci = 50, traits = 4), mutation(0.0005, 0.2), env_sd = 0.8,
                      recombination = genetic_map(chromosomes = 2, length = 0.5))
pop <- evolve(model, stabilizing(2000, omega = 10), n = 1000, record_every = 50, seed = 1)
g_decomposition(pop)           # G, genic and ld for the current population
g_matrices(pop, "ld")          # the LD part at every recorded generation
```

Pleiotropy shows up in the genic part (one allele affecting several traits),
while correlations built by selection between loci show up in the LD part.

For any architecture, including `multilinear()`, the split can be estimated
by shuffling alleles among individuals at each locus, which removes
associations between loci but keeps allele frequencies (as in do O & Whitlock
2023):

```r
g_decomposition(pop, method = "shuffle", shuffles = 20)
```

### Mutational covariance

`mutational_covariance()` applies single mutations to genomes drawn from the
population and measures their effects on genotypic values. It returns the
mean effect, the covariance of the effects of one mutation, and the
mutational (co)variance input per generation, `M`. With epistasis or a
binary B matrix, mutational effects depend on the genetic background, so `M`
itself evolves:

```r
mutational_covariance(pop, n = 10000)
```

### The model

* **Genomes**: each allele carries an effect vector in trait space; the effect
  of a locus, `Y_l`, is the sum of its two alleles. Effects are measured from
  the reference genotype (all effects zero), where new populations start
  unless `init_sd > 0`.
* **Architectures**:
  * `additive()`: `x = Σ Y_l`.
  * `multilinear()`: `x = S + Σ_{i<j} ε(Y_i, Y_j)` with `S = Σ Y_l` and a
    `traits × traits × traits` tensor `ε` shared by all pairs of loci
    (Hansen & Wagner 2001).
  * `binary_b()`: the model of Melo & Marroig (2015). Each allele has a
    scalar effect `y` and a column of a binary matrix `B` saying which traits
    it affects, `x = B y`. Used with `binary_b_mutation()`, which mutates `y`
    and the entries of `B` separately (`b_mutation = "random"` reproduces the
    B mutation of the published simulations, which changed entries at half
    the nominal rate).

  `local_effects()` gives each locus's effect in the background of the rest
  of the genome.
* **Phenotypes**: `z = x + e`, `e ~ N(0, env_sd²)`.
* **Recombination**: loci are unlinked by default. `genetic_map(chromosomes,
  length)` places loci evenly along chromosomes (Haldane map function), and
  `linkage(rates)` sets the recombination fraction between each pair of
  neighbouring loci; pass either as `evomod_model(..., recombination = )`.
* **Each generation**: the optimum moves (directional regimes), fitness
  `w = exp(-(z - θ)' Ω⁻¹ (z - θ) / 2)` is computed, every individual mutates
  (each allele with probability `rate`, adding a `N(0, cov)` vector), and `n`
  offspring are produced from parents drawn with probability proportional to
  fitness, with free recombination.
* **Random numbers** come from R's generator: `set.seed()` or `seed` makes
  runs reproducible, and continuing a population gives the same result as a
  single run.

## Extending evomod2

The engine is a set of C++ components declared in `src/evomod.h`:

| Interface | Role | Implementations |
|---|---|---|
| `Architecture` | genome → genotypic value, local effects | `Additive`, `Multilinear`, `BinaryB` |
| `Mutation` | changes allelic effects | `GaussianMutation`, `BinaryBMutation` |
| `Fitness` | phenotypes → log fitness | `GaussianFitness` |
| `Recombination` | parents → offspring genome | `FreeRecombination`, `LinkageMap` |

To add, for example, a new architecture:

1. Derive a class from `Architecture` in `src/evomod.h` and implement
   `genotypic_value()` and `local_effects()` in `src/components.cpp`.
2. Add a case to `make_architecture()` that builds it from the list the R
   constructor returns.
3. Write the R constructor in `R/architecture.R`, returning
   `structure(list(type = "<name>", loci, traits, state_rows, mutation_type, ...), class = "evomod_architecture")`.
   `state_rows` is the number of rows of a genome (one column per allele):
   `traits` for an effect vector per allele, `traits + 1` for `binary_b()`.
   `mutation_type` is the mutation model the architecture works with; if the
   genome needs a special starting state, add it to `init_genomes()`.
4. Add tests in `tests/testthat/` (compare with an explicit, slow R version of
   the same map), and, for new models, prototype and cross-check against
   EvoMod.jl.

## Validation

The scripts in `validation/` run matching simulations in evomod2 and in an
independent implementation, and compare replicate means of G matrices and
trait means:

```sh
# additive and multilinear models against EvoMod.jl
Rscript validation/cross_check.R path/to/evomod_julia 20
# binary B model against the original C code (lem-usp/evomod, c-gsl; needs GSL)
Rscript validation/cross_check_c.R path/to/evomod/c-gsl 100 8
```

The C code moves the optimum after each generation and skips the first
directional step, so under directional selection its trait means lag evomod2's
by about two generations of optimum movement.

## References

Hansen, T. F., & Wagner, G. P. (2001). Modeling genetic architecture: a
multilinear theory of gene interaction. *Theoretical Population Biology*,
59(1), 61-86.

Melo, D., & Marroig, G. (2015). Directional selection can drive the evolution
of modularity in complex traits. *PNAS*, 112(2), 470-475.
