#' Recorded population statistics
#'
#' Statistics recorded during [evolve()], every `record_every` generations.
#'
#' * `trait_means()`: mean genotypic value (`mean_x`), mean phenotype
#'   (`mean_z`) and optimum (`theta`) of each trait, in long format.
#' * `variances()`: additive genetic (`G`) and phenotypic (`P`) variances
#'   and heritability (`h2`) of each trait, in long format.
#' * `g_matrices()`, `p_matrices()`: lists of covariance matrices of genotypic
#'   values (G) and phenotypes (P), named by generation. For architectures
#'   where genotypic values are sums of allele effects ([additive()],
#'   [binary_b()]) G can be split into its genic part, the G expected
#'   without linkage disequilibrium (twice the sum over loci of the covariance
#'   of allele effects), and the rest, due to linkage disequilibrium between
#'   loci and departures from Hardy-Weinberg (`component = "ld"`); see also
#'   [g_decomposition()].
#'
#' @param pop an `evomod_population`.
#' @param component `"total"` (G), `"genic"` or `"ld"`.
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
  a <- if (which == "G_ld") pop$records$G - pop$records$G_genic else pop$records[[which]]
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
g_matrices <- function(pop, component = c("total", "genic", "ld")) {
  component <- match.arg(component)
  if (component != "total" && is.null(pop$records$G_genic))
    stop("the genic part of G is only defined for architectures where genotypic values ",
         "are sums of allele effects (additive, binary_b)", call. = FALSE)
  matrices(pop, c(total = "G", genic = "G_genic", ld = "G_ld")[[component]])
}

#' @rdname statistics
#' @export
p_matrices <- function(pop) matrices(pop, "P")

#' Genomes and genotypic values
#'
#' * `genomes()`: genomes of the current population, an array
#'   `rows x alleles x individuals` (alleles `2l - 1` and `2l` are locus `l`).
#'   Rows are the allelic effect vector (one per trait), or, for [binary_b()],
#'   `y` followed by the column of `B`.
#' * `genotypic_values()`: genotypic values of genomes under an architecture
#'   (traits x individuals).
#' * `local_effects()`: effect of each locus of one genome in the background of
#'   the rest of the genome (traits x loci): the change in genotypic value when
#'   the locus is set to the reference genotype. Equal to the locus effects
#'   for the additive architecture.
#'
#' @param pop an `evomod_population`.
#' @param x an architecture, model or population.
#' @param alleles a genome (rows x alleles) or an array of genomes
#'   (rows x alleles x individuals), see `genomes()`.
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
  check_genome_rows(alleles, arch)
  genotypic_values_cpp(alleles, arch$traits, arch)
}

#' @rdname genomes
#' @export
local_effects <- function(x, alleles = genomes(x)[, , individual], individual = 1) {
  arch <- as_architecture(x)
  alleles <- as.matrix(alleles)
  check_genome_rows(alleles, arch)
  local_effects_cpp(alleles, arch)
}

check_genome_rows <- function(alleles, arch) {
  if (dim(alleles)[1] != arch$state_rows || dim(alleles)[2] != 2 * arch$loci)
    stop("genomes of this architecture are ", arch$state_rows, " x ", 2 * arch$loci,
         " (rows x alleles)", call. = FALSE)
}

#' Genic and linkage disequilibrium parts of G
#'
#' Splits the G matrix of the current population into its genic part, the G
#' expected without linkage disequilibrium or departures from Hardy-Weinberg
#' (twice the sum over loci of the covariance of allele effects, pooling both
#' alleles of every individual), and the remainder (`ld`), due to linkage
#' disequilibrium between loci and departures from Hardy-Weinberg. Defined for
#' architectures where genotypic values are sums of allele effects
#' ([additive()], [binary_b()]). All matrices use the `n - 1` denominator of
#' [stats::cov()], so `G = genic + ld`.
#'
#' @param pop an `evomod_population`.
#' @return A list with matrices `G`, `genic` and `ld`.
#' @examples
#' model <- evomod_model(additive(loci = 20, traits = 2),
#'                       mutation(rate = 0.01, cov = 0.1), init_sd = 0.1)
#' pop <- evolve(model, stabilizing(20, omega = 10), n = 100, seed = 1)
#' g_decomposition(pop)
#' @export
g_decomposition <- function(pop) {
  arch <- pop$model$architecture
  genic <- genic_covariance_cpp(genomes(pop), arch)
  if (!length(genic))
    stop("the genic part of G is only defined for architectures where genotypic values ",
         "are sums of allele effects (additive, binary_b)", call. = FALSE)
  G <- stats::cov(t(pop$x))
  nm <- trait_names(pop)
  dimnames(genic) <- dimnames(G) <- list(nm, nm)
  list(G = G, genic = genic, ld = G - genic)
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
