#' evomod2: individual-based simulations of genetic architecture evolution
#'
#' Build a model from an architecture ([additive()], [multilinear()]) and a
#' mutation model ([mutation()]) with [evomod_model()], then run it through
#' selection regimes ([drift()], [stabilizing()], [directional()]) with
#' [evolve()]. Recorded statistics: [trait_means()], [variances()],
#' [g_matrices()].
#'
#' New components (architectures, mutation models, fitness functions) are
#' C++ classes in `src/evomod.h`; see the "Extending evomod2" section of the
#' README.
#'
#' @useDynLib evomod2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"
