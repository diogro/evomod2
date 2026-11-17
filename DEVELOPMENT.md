# DEVELOPMENT.md

Guidance for AI coding agents (and people) extending or modifying
**evomod2**. For using the package, see [AGENTS.md](AGENTS.md); the
conventions listed there (standard deviations, per-allele rates, genome
layout, optimum carry-over) apply to new code too.

## Layout

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

## The generation loop

`evolve_cpp()` in `src/simulate.cpp` runs the regimes. In every generation, in
this order:

1. the optimum moves (`directional()`: `theta += delta_theta`;
   `fluctuating()`: a jump after every `every` generations);
2. fitness is computed from the current phenotypes (all ones under drift);
3. reproduction:
   - Wright-Fisher: every genome mutates (`Mutation::mutate`), `n` sires and
     `n` dames are drawn with probability proportional to fitness, and
     `Recombination::offspring` builds the offspring into a separate buffer,
     which is then swapped in (parents are never overwritten);
   - Moran: `n` events; each replaces a random individual by an offspring of
     two parents drawn from a Fenwick tree of fitnesses, the offspring then
     mutates, and only its fitness is updated;
4. genotypic values (`Architecture::genotypic_value`) and phenotypes
   (`x + Normal(0, env_sd)`) are computed;
5. statistics are recorded when `generation % record_every == 0`.

The optimum carries over between regimes; a regime resets it only when the
user gives `theta` (`reset_theta`).

## Interfaces in detail

- `Architecture`: `genotypic_value(genome)` (traits vector),
  `local_effects(genome)` (traits x loci: the change in genotypic value when
  a locus is set to the reference genotype, all zeros), and optionally
  `allele_effects(genome)` (traits x alleles) when genotypic values are sums
  of allele effects; return an empty matrix otherwise. `allele_effects`
  enables the exact `g_decomposition()` and the recorded genic part of G.
- `Mutation`: `mutate(genome)` (one generation of mutation, using binomial
  numbers of events), `mutate_one(genome)` (exactly one event, at a site
  chosen in proportion to the rates; used by `mutational_covariance()`), and
  `expected_mutations(genome)` (events per genome per generation). Support
  `persistence` (`new = persistence * old + draw`) and keep the `+=` path
  when it is 1.
- `Recombination`: `offspring(sire, dame, kid)`; column `2l` of the kid gets
  one of the sire's alleles at locus `l` (0-based), column `2l + 1` one of
  the dam's.
- `Fitness`: `log_fitness(Z)` for a traits x individuals matrix.

## Random numbers

Use only R's generator from C++ (`R::unif_rand()`, `R::norm_rand()`,
`R::rbinom()`, `rnorm_vec()`), never `<random>` or Armadillo's generator, so
that `set.seed()` and `seed =` reproduce runs and a population continued
over several `evolve()` calls matches a single call. Draw the random numbers
in the same order on every code path you do not mean to change.

## Adding a regime or a statistic

- A regime is an R constructor in `R/regimes.R` returning a `new_regime()`
  list; `resolve_regimes()` turns it into the list `evolve_cpp()` reads
  (`generations`, `selective`, `reset_theta`, `theta`, `delta_theta`,
  `omega`, and regime-specific fields). Add the per-generation behaviour to
  the loop in `evolve_cpp()`.
- A statistic recorded during the run goes into `Records` in
  `src/simulate.cpp` (and into `append_records()` in `R/evolve.R` so that
  continued runs keep it), with an accessor in `R/accessors.R`. Statistics
  of the current population that are expensive or optional (shuffles,
  mutational covariance, additive covariance) are better computed on demand
  from `genomes(pop)` in R or through a C++ helper.

## Adding a component

Follow the README section "Extending evomod2": derive a class in
`src/evomod.h`, implement it in `src/components.cpp`, add a case to the
matching `make_*()` factory, write the R constructor returning a list with a
`type` (architectures also need `loci`, `traits`, `state_rows`,
`mutation_type`), and test it against an explicit, slow R version of the same
computation. New models are usually prototyped and cross-checked in EvoMod.jl
first.

## Build, test, document

```sh
Rscript -e 'Rcpp::compileAttributes(); roxygen2::roxygenise()'   # after changing exports or docs
R CMD INSTALL --preclean .
Rscript -e 'devtools::test()'
R CMD build . && R CMD check --no-manual evomod2_*.tar.gz
Rscript pkgdown/build_site.R                                     # website into docs/
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

## Validation, CI and the website

- `validation/cross_check.R` (against EvoMod.jl, additive and multilinear)
  and `validation/cross_check_c.R` (against the original C code, binary B)
  compare replicate means of G and trait means; run them after changes to the
  generation loop, mutation or the architectures they cover.
- GitHub Actions runs `R CMD check` on macOS, Ubuntu and Windows for every
  push to `main` and every pull request.
- The website is built by pkgdown into `docs/` and served by GitHub Pages
  from `main`; rebuild it (`Rscript pkgdown/build_site.R`) when the README,
  documentation or vignettes change. The home page is the README. Add new
  vignettes to `_pkgdown.yml`. `pkgdown/build_site.R` and
  `pkgdown/extra.scss` are shared with the sibling packages' sites; keep the
  copies identical.

