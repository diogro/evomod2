# AGENTS.md

Guidance for AI coding agents helping people use or develop **evomod2**, an R
package (C++ core via Rcpp/RcppArmadillo) for individual-based simulations of
multivariate trait evolution: genetic architectures, mutation, recombination,
selection regimes and G/M matrix statistics.

Documentation: <https://diogro.github.io/evomod2/>. Start from the vignettes;
they are the best source of working code.

## Using the package

### Install

```r
remotes::install_github("diogro/evomod2", build_vignettes = TRUE)
```

A C++ compiler is required (Xcode command line tools on macOS, Rtools on
Windows).

### The workflow

Every simulation is a model, a list of regimes, and `evolve()`:

```r
library(evomod2)
model <- evomod_model(
  architecture = additive(loci = 50, traits = 4),
  mutation = mutation(rate = 0.0005, cov = 0.2),   # cov: a matrix, or a number = sd per trait
  env_sd = 0.8,
  init_sd = 0.1,
  recombination = free_recombination(),           # default
  life_cycle = wright_fisher())                   # default
pop <- evolve(model, n = 1000, record_every = 50, seed = 1,
              regimes = list(drift(500),
                             stabilizing(2000, omega = 10),
                             directional(2000, delta_theta = c(0.01, 0.01, 0, 0), omega = 10)))
pop <- evolve(pop, stabilizing(500, omega = 10))  # continue the same population
```

Read results with `trait_means()`, `variances()`, `g_matrices()`,
`p_matrices()`, `mean_fitness()`, `plot(pop)`. Analyse with `avg_corr()`,
`autonomy()`, `mean_autonomy()`, `eigen_variance()`, `g_decomposition()`,
`additive_covariance()`, `mutational_covariance()`, `local_effects()`.

### Components (pick one of each)

| Component | Options |
|---|---|
| Architecture | `additive()`, `binary_b()`, `multilinear()` (shared tensor from `random_epsilon()`, or one tensor per pair of loci from `pairwise_epsilon()`), `contribution_matrix()` |
| Mutation | `mutation()` / `hoc_mutation()` for additive and multilinear; `binary_b_mutation()` for `binary_b()`; `contribution_mutation()` for `contribution_matrix()` (the model checks the match) |
| Recombination | `free_recombination()`, `genetic_map(chromosomes, length)`, `linkage(rates)` |
| Life cycle | `wright_fisher()` (discrete generations), `moran()` (overlapping) |
| Regimes | `drift()`, `stabilizing()`, `directional()`, `fluctuating()` |

### Conventions that are easy to get wrong

- **Standard deviations, not variances.** `env_sd`, `init_sd`, the `sd`
  arguments and a scalar mutational `cov` are standard deviations; a matrix
  `cov` is a covariance matrix. `omega` is a covariance matrix (a scalar
  means `omega * I`). Fitness is `exp(-(z - theta)' omega^-1 (z - theta) / 2)`.
- **Rates are per allele (or per entry of B) per generation.** A diploid
  locus has two alleles, so the expected number of mutations per genome is
  `rate * 2 * loci` (see `mutational_covariance(pop)$expected_mutations`).
- **Genomes** are arrays `rows x alleles x individuals`
  (`genomes(pop)`); alleles `2l - 1` and `2l` are the paternal and maternal
  copies of locus `l`. Rows are the effect vector (one per trait), except for
  `binary_b()` (`y`, then the column of `B`) and `contribution_matrix()` (one
  scalar).
- **The optimum carries over** between regimes and between `evolve()` calls;
  a regime resets it only when given `theta`. `directional()` moves it before
  selection in every generation; `fluctuating()` jumps after every `every`
  generations.
- **Records** are taken every `record_every` generations (and at the start of
  a new population); the first record is generation 0. Moments are recorded
  after the generation's reproduction.
- **`persistence`** in the mutation functions: 1 (default) is a random walk,
  0 the house of cards model.
- **`g_decomposition()`** is exact only for architectures whose genotypic
  values are sums of allele effects (`additive()`, `binary_b()`); for the
  others use `method = "shuffle"`.
- **`additive_covariance()`** needs more individuals than loci x traits.
- **Moran drift** is about twice as strong per generation as Wright-Fisher at
  the same `n` (`N_e` is about `n / 2`).
- **Binary B mutation:** `b_mutation = "flip"` (default) always changes a
  mutated entry; `"random"` reproduces the original C code of Melo & Marroig
  (2015), where entries changed at half the nominal rate.

### Reproducibility

Random numbers come from R's generator: `seed = ` (or `set.seed()`) makes a run
reproducible, and continuing a population in several `evolve()` calls gives
the same result as one call. When running replicates in parallel, give each
its own seed.

### Running many replicates

The examples in `system.file("examples", package = "evomod2")` run replicated
experiments from the command line, one replicate per process
(`Rscript <script> <mode> ... <outdir>`, then a `summarise` mode). Use that
pattern for large experiments: one single-threaded process per replicate,
started with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1` so
that multithreaded BLAS does not oversubscribe the machine. Typical costs: a
few seconds for 1000 individuals with 20 loci and 4 traits over 5000
generations; about 0.25 s per generation for 5000 individuals with 500 loci
and 10 traits (2 GB of memory).

### Published results reproduced in the vignettes

- `vignette("melo-marroig-2015")`: directional selection and modularity.
- `vignette("fluctuating-environments")`: do O & Whitlock (2023) across
  architectures.
- `vignette("jones-burger-arnold-2014")`: epistasis and the evolution of M.

Each states the parameters used and where evomod2 differs from the original
simulations. Point users to them before re-deriving those setups.

## Developing the package

### Layout

- `src/evomod.h`: the C++ component interfaces (`Architecture`, `Mutation`,
  `Fitness`, `Recombination`) and their implementations; `src/components.cpp`
  builds components from the R lists (`make_*()` factories);
  `src/simulate.cpp` has the generation loop (`evolve_cpp`) and C++ helpers
  exported to R.
- `R/`: constructors (`architecture.R`, `model.R`, `regimes.R`,
  `recombination.R`), `evolve.R`, accessors and analysis (`accessors.R`,
  `statistics.R`), `plot.R`.
- `tests/testthat/`: unit tests; `validation/`: cross-checks against
  independent implementations (EvoMod.jl, the original C code), not run in CI.
- `inst/examples/`: command-line experiments; `inst/extdata/`: stored
  summaries used by the results vignettes; `vignettes/`; `_pkgdown.yml` and
  `docs/` (the website).

### Adding a component

Follow the README section "Extending evomod2": derive a class in
`src/evomod.h`, implement it in `src/components.cpp`, add a case to the
matching `make_*()` factory, write the R constructor returning a list with a
`type` (architectures also need `loci`, `traits`, `state_rows`,
`mutation_type`), and test it against an explicit, slow R version of the same
computation. New models are usually prototyped and cross-checked in EvoMod.jl
first.

### Build, test, document

```sh
Rscript -e 'Rcpp::compileAttributes(); roxygen2::roxygenise()'   # after changing exports or docs
R CMD INSTALL --preclean .
Rscript -e 'devtools::test()'
R CMD build . && R CMD check --no-manual evomod2_*.tar.gz
Rscript -e 'pkgdown::build_site()'                               # website into docs/
```

- **Time code only from clean builds.** `devtools::test()` and `load_all()`
  compile without optimisation and leave the objects in `src/`; a later
  `R CMD INSTALL .` without `--preclean` reuses them and runs about three
  times slower. Remove `src/*.o` and `src/*.so` after testing.
- Keep existing behaviour bit-for-bit when adding options with defaults
  (seeded runs before and after the change must be identical); compilers may
  fuse multiply-adds, so keep the old arithmetic on the default path.
- Statistical tests: use enough replicates and tolerances that they do not
  fail by chance, and prefer exact expectations (identities, analytical
  equilibria) where they exist.
- Documentation is roxygen2 (markdown); do not edit `man/` or `NAMESPACE`
  by hand.
- Vignettes that run simulations must build in under a minute; heavy results
  go in `inst/extdata/` with the script that produces them in
  `inst/examples/`.
