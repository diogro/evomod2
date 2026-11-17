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

### Documentation

* `vignette("evomod2")`: getting started, models, regimes, recorded
  statistics, G decomposition.
* `vignette("architectures")`: binary B, additive and multilinear
  architectures; epistasis; linkage and the genic/LD parts of G.
* `vignette("mutation-models")`: random walk and house of cards mutation
  (`persistence`, `hoc_mutation()`): neutral equilibria, stabilizing
  selection and a moving optimum.
* `vignette("melo-marroig-2015")`: reproducing the simulations of Melo &
  Marroig (2015), with results of a full run.
* `vignette("jones-burger-arnold-2014")`: reproducing Jones, Bürger &
  Arnold (2014) on the evolution of the mutational matrix M under epistasis,
  with results of a full run.
* `vignette("fluctuating-environments")`: the claims of do O & Whitlock
  (2023) about fluctuating optima, tested on four genetic architectures,
  with results of a full run.
* Command-line examples in `system.file("examples", package = "evomod2")`:
  `corridor_selection.R`, `corridor_linkage.R`, `melo_marroig_2015.R`,
  `fluctuating_environments.R`, `jones_burger_arnold_2014.R`.

Install with `remotes::install_github("diogro/evomod2", build_vignettes = TRUE)`
to get the vignettes.

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
  * `contribution_matrix()`: the model of do O & Whitlock (2023), `x = C a`
    with weight genes `a` and contribution genes `C` as separate loci, used
    with `contribution_mutation()`.
  * `binary_b()`: the model of Melo & Marroig (2015). Each allele has a
    scalar effect `y` and a column of a binary matrix `B` saying which traits
    it affects, `x = B y`. Used with `binary_b_mutation()`, which mutates `y`
    and the entries of `B` separately (`b_mutation = "random"` reproduces the
    B mutation of the published simulations, which changed entries at half
    the nominal rate).

  `local_effects()` gives each locus's effect in the background of the rest
  of the genome.
* **Phenotypes**: `z = x + e`, `e ~ N(0, env_sd²)`.
* **Life cycle**: discrete generations (`wright_fisher()`, the default) or
  overlapping generations (`moran()`: a generation is `n` birth-death events,
  each replacing a random individual by the offspring of two parents drawn
  in proportion to fitness), or the life cycle of Jones et al.
  (`viability_selection()`: separate sexes, females mating with several
  males, viability selection on zygotes); `evomod_model(..., life_cycle = moran())`.
* **Mutation**: a mutation replaces an allele's effect `e` by
  `persistence * e + draw`; `persistence = 1` (default) is a random walk,
  `0` the house of cards model (`hoc_mutation()`), for every architecture.
* **Recombination**: loci are unlinked by default. `genetic_map(chromosomes,
  length)` places loci evenly along chromosomes (Haldane map function), and
  `linkage(rates)` sets the recombination fraction between each pair of
  neighbouring loci; pass either as `evomod_model(..., recombination = )`.
* **Regimes**: `drift()`, `stabilizing()`, `directional()` (the optimum moves
  steadily) and `fluctuating()` (the optimum jumps up or down every few
  generations, all traits together, independently, or in blocks; do O &
  Whitlock 2023).
* **Each generation**: the optimum moves (directional and fluctuating regimes), fitness
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

### Models implemented

- **Binary B matrix:** Wagner, G. P. (1989). Multivariate mutation-selection
  balance with constrained pleiotropic effects. *Genetics*, 122(1), 223-234.
  [doi:10.1093/genetics/122.1.223](https://doi.org/10.1093/genetics/122.1.223);
  Melo, D., & Marroig, G. (2015). Directional selection can drive the
  evolution of modularity in complex traits. *PNAS*, 112(2), 470-475.
  [doi:10.1073/pnas.1322632112](https://doi.org/10.1073/pnas.1322632112)
- **Multilinear epistasis:** Hansen, T. F., & Wagner, G. P. (2001). Modeling
  genetic architecture: a multilinear theory of gene interaction. *Theoretical
  Population Biology*, 59(1), 61-86.
  [doi:10.1006/tpbi.2000.1508](https://doi.org/10.1006/tpbi.2000.1508);
  Hermisson, J., Hansen, T. F., & Wagner, G. P. (2003). Epistasis in polygenic
  traits and the evolution of genetic architecture under stabilizing
  selection. *The American Naturalist*, 161(5), 708-734.
  [doi:10.1086/374204](https://doi.org/10.1086/374204)
- **Additive pleiotropic mutation, random epistasis between pairs of loci, and
  the `viability_selection()` life cycle:** Jones, A. G., Arnold, S. J., &
  Bürger, R. (2003). Stability of the G-matrix in a population experiencing
  pleiotropic mutation, stabilizing selection, and genetic drift. *Evolution*,
  57(8), 1747-1760. [doi:10.1554/02-631](https://doi.org/10.1554/02-631);
  Jones, A. G., Arnold, S. J., & Bürger, R. (2004). Evolution and stability of
  the G-matrix on a landscape with a moving optimum. *Evolution*, 58(8),
  1639-1654. [doi:10.1554/03-651](https://doi.org/10.1554/03-651);
  Jones, A. G., Bürger, R., & Arnold, S. J. (2014). Epistasis and natural
  selection shape the mutational architecture of complex traits. *Nature
  Communications*, 5, 3709.
  [doi:10.1038/ncomms4709](https://doi.org/10.1038/ncomms4709)
- **Contribution matrix and fluctuating optima:** do O, I., & Whitlock, M. C.
  (2023). The evolution of genetic covariance and modularity as a result of
  multigenerational environmental fluctuation. *Evolution Letters*, 7(6),
  457-466. [doi:10.1093/evlett/qrad048](https://doi.org/10.1093/evlett/qrad048)
- **House of cards mutation:** Kingman, J. F. C. (1978). A simple model for
  the balance between selection and mutation. *Journal of Applied
  Probability*, 15(1), 1-12.
  [doi:10.2307/3213231](https://doi.org/10.2307/3213231);
  Turelli, M. (1984). Heritable genetic variation via mutation-selection
  balance: Lerch's zeta meets the abdominal bristle. *Theoretical Population
  Biology*, 25(2), 138-193.
  [doi:10.1016/0040-5809(84)90017-0](https://doi.org/10.1016/0040-5809(84)90017-0);
  Zeng, Z.-B., & Cockerham, C. C. (1993). Mutation models and quantitative
  genetic variation. *Genetics*, 133(3), 729-736.
  [doi:10.1093/genetics/133.3.729](https://doi.org/10.1093/genetics/133.3.729)
- **Genetic map:** Haldane, J. B. S. (1919). The combination of linkage
  values and the calculation of distances between the loci of linked factors.
  *Journal of Genetics*, 8(4), 299-309.

### Results reproduced in the vignettes

- Melo & Marroig (2015), `vignette("melo-marroig-2015")`
- Jones, Bürger & Arnold (2014), `vignette("jones-burger-arnold-2014")`
- do O & Whitlock (2023), `vignette("fluctuating-environments")`

### Theory and statistics

- **Response to selection and drift of multivariate traits:** Lande, R.
  (1976). Natural selection and random genetic drift in phenotypic evolution.
  *Evolution*, 30(2), 314-334.
  [doi:10.1111/j.1558-5646.1976.tb00911.x](https://doi.org/10.1111/j.1558-5646.1976.tb00911.x);
  Lande, R. (1979). Quantitative genetic analysis of multivariate evolution,
  applied to brain: body size allometry. *Evolution*, 33(1), 402-416.
  [doi:10.2307/2407630](https://doi.org/10.2307/2407630)
- **Linkage disequilibrium under selection (genic and LD parts of G):**
  Bulmer, M. G. (1971). The effect of selection on genetic variability.
  *The American Naturalist*, 105(943), 201-211.
  [doi:10.1086/282718](https://doi.org/10.1086/282718)
- **Evolvability and autonomy:** Hansen, T. F., & Houle, D. (2008). Measuring
  and comparing evolvability and constraint in multivariate characters.
  *Journal of Evolutionary Biology*, 21(5), 1201-1219.
  [doi:10.1111/j.1420-9101.2008.01573.x](https://doi.org/10.1111/j.1420-9101.2008.01573.x)
- **Modularity:** Melo, D., Porto, A., Cheverud, J. M., & Marroig, G. (2016).
  Modularity: genes, development, and evolution. *Annual Review of Ecology,
  Evolution, and Systematics*, 47, 463-486.
  [doi:10.1146/annurev-ecolsys-121415-032409](https://doi.org/10.1146/annurev-ecolsys-121415-032409)
