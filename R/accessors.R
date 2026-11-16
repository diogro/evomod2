#' Recorded population statistics
#'
#' Statistics recorded during [evolve()], every `record_every` generations.
#'
#' * `trait_means()`: mean genotypic value (`mean_x`), mean phenotype
#'   (`mean_z`) and optimum (`theta`) of each trait, in long format.
#' * `variances()`: additive genetic (`G`) and phenotypic (`P`) variances
#'   and heritability (`h2`) of each trait, in long format.
#' * `g_matrices()`, `p_matrices()`: lists of covariance matrices of genotypic
#'   values (G) and phenotypes (P), named by generation.
#'
#' @param pop an `evomod_population`.
#' @return A data frame or a named list of matrices.
#' @name statistics
NULL

trait_names <- function(pop) paste0("trait", seq_len(pop$model$traits))

#' @rdname statistics
#' @export
trait_means <- function(pop) {
  r <- pop$records
  p <- pop$model$traits
  data.frame(generation = rep(r$generation, each = p),
             trait = rep(trait_names(pop), length(r$generation)),
             mean_x = as.vector(r$mean_x), mean_z = as.vector(r$mean_z),
             theta = as.vector(r$theta))
}

#' @rdname statistics
#' @export
variances <- function(pop) {
  r <- pop$records
  p <- pop$model$traits
  g <- apply(r$G, 3, diag)
  pp <- apply(r$P, 3, diag)
  data.frame(generation = rep(r$generation, each = p),
             trait = rep(trait_names(pop), length(r$generation)),
             G = as.vector(g), P = as.vector(pp), h2 = as.vector(g / pp))
}

matrices <- function(pop, which) {
  a <- pop$records[[which]]
  out <- lapply(seq_len(dim(a)[3]), function(i) {
    m <- a[, , i, drop = TRUE]
    m <- matrix(m, pop$model$traits, pop$model$traits)
    dimnames(m) <- list(trait_names(pop), trait_names(pop))
    m
  })
  names(out) <- pop$records$generation
  out
}

#' @rdname statistics
#' @export
g_matrices <- function(pop) matrices(pop, "G")

#' @rdname statistics
#' @export
p_matrices <- function(pop) matrices(pop, "P")

#' Genomes and genotypic values
#'
#' * `genomes()`: allelic effects of the current population, an array
#'   `traits x alleles x individuals` (alleles `2l - 1` and `2l` are locus `l`).
#' * `genotypic_values()`: genotypic values of genomes under an architecture
#'   (traits x individuals).
#' * `local_effects()`: effect of each locus of one genome in the background of
#'   the rest of the genome (traits x loci): the change in genotypic value when
#'   the locus is set to the reference genotype. Equal to the locus effects
#'   for the additive architecture.
#'
#' @param pop an `evomod_population`.
#' @param x an architecture, model or population.
#' @param alleles a genome (traits x alleles) or an array of genomes
#'   (traits x alleles x individuals).
#' @param individual index of an individual of `pop`.
#' @name genomes
NULL

#' @rdname genomes
#' @export
genomes <- function(pop) pop$alleles

as_architecture <- function(x) {
  if (inherits(x, "evomod_population")) x <- x$model
  if (inherits(x, "evomod_model")) x <- x$architecture
  if (!inherits(x, "evomod_architecture")) stop("`x` must be an architecture, model or population")
  x
}

#' @rdname genomes
#' @export
genotypic_values <- function(x, alleles = genomes(x)) {
  arch <- as_architecture(x)
  if (length(dim(alleles)) == 2) alleles <- array(alleles, c(dim(alleles), 1))
  genotypic_values_cpp(alleles, arch)
}

#' @rdname genomes
#' @export
local_effects <- function(x, alleles = genomes(x)[, , individual], individual = 1) {
  local_effects_cpp(as.matrix(alleles), as_architecture(x))
}

#' Average correlations within and between modules
#'
#' @param m a correlation (or covariance) matrix.
#' @param modules list of trait indices of each module; by default the first
#'   and second half of the traits.
#' @return A list with the mean correlation within each module (`within`), the
#'   mean correlation between traits in different modules (`between`), and
#'   their ratio (`avg_ratio`).
#' @examples
#' m <- matrix(c(1, .8, .1, .1, .8, 1, .1, .1, .1, .1, 1, .6, .1, .1, .6, 1), 4)
#' avg_corr(m)
#' @export
avg_corr <- function(m, modules = NULL) {
  m <- stats::cov2cor(as.matrix(m))
  p <- nrow(m)
  if (is.null(modules)) modules <- list(seq_len(p %/% 2), (p %/% 2 + 1):p)
  within <- vapply(modules, function(mod) mean(m[mod, mod][lower.tri(diag(length(mod)))]), 1)
  same <- matrix(FALSE, p, p)
  for (mod in modules) same[mod, mod] <- TRUE
  between <- mean(m[!same])
  list(within = within, between = between, avg_ratio = mean(within) / between)
}
