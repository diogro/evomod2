# evomod2: individual-based simulations of genetic architecture evolution

Build a model from an architecture
([`additive()`](https://diogro.github.io/evomod2/reference/architectures.md),
[`multilinear()`](https://diogro.github.io/evomod2/reference/architectures.md),
[`binary_b()`](https://diogro.github.io/evomod2/reference/architectures.md))
and a mutation model
([`mutation()`](https://diogro.github.io/evomod2/reference/mutation.md),
[`binary_b_mutation()`](https://diogro.github.io/evomod2/reference/mutation.md))
with
[`evomod_model()`](https://diogro.github.io/evomod2/reference/evomod_model.md),
then run it through selection regimes
([`drift()`](https://diogro.github.io/evomod2/reference/regimes.md),
[`stabilizing()`](https://diogro.github.io/evomod2/reference/regimes.md),
[`directional()`](https://diogro.github.io/evomod2/reference/regimes.md))
with [`evolve()`](https://diogro.github.io/evomod2/reference/evolve.md).
Recorded statistics:
[`trait_means()`](https://diogro.github.io/evomod2/reference/statistics.md),
[`variances()`](https://diogro.github.io/evomod2/reference/statistics.md),
[`g_matrices()`](https://diogro.github.io/evomod2/reference/statistics.md),
[`g_decomposition()`](https://diogro.github.io/evomod2/reference/g_decomposition.md).
Loci are unlinked unless a
[`genetic_map()`](https://diogro.github.io/evomod2/reference/recombination.md)
is given to
[`evomod_model()`](https://diogro.github.io/evomod2/reference/evomod_model.md).

## Details

New components (architectures, mutation models, fitness functions) are
C++ classes in `src/evomod.h`; see the "Extending evomod2" section of
the README.

## See also

Useful links:

- <https://diogro.github.io/evomod2/>

- <https://github.com/diogro/evomod2>

- Report bugs at <https://github.com/diogro/evomod2/issues>

## Author

**Maintainer**: Diogo Melo <diogro@gmail.com>

Authors:

- Diogo Melo <diogro@gmail.com>
