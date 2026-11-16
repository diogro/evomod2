#' Recorded population statistics
#'
#' Statistics recorded during [evolve()], every `record_every` generations.
#'
#' * `trait_means()`: mean genotypic value (`mean_x`), mean phenotype
#'   (`mean_z`) and optimum (`theta`) of each trait, in long format.
#' * `mean_fitness()`: mean fitness of the population on the fitness surface
#'   of the current regime (`NA` under drift).
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
mean_fitness <- function(pop) {
  data.frame(generation = pop$records$generation, mean_fitness = pop$records$mean_fitness)
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
local_effects <- function(x, alleles = genome_of(x, individual), individual = 1) {
  arch <- as_architecture(x)
  alleles <- as.matrix(alleles)
  check_genome_rows(alleles, arch)
  local_effects_cpp(alleles, arch)
}

genome_of <- function(pop, i) {
  g <- genomes(pop)
  matrix(g[, , i], nrow = dim(g)[1])
}

check_genome_rows <- function(alleles, arch) {
  if (dim(alleles)[1] != arch$state_rows || dim(alleles)[2] != 2 * arch$loci)
    stop("genomes of this architecture are ", arch$state_rows, " x ", 2 * arch$loci,
         " (rows x alleles)", call. = FALSE)
}

#' Pleiotropic and linkage disequilibrium parts of G
#'
#' Splits the G matrix of the current population into the G expected without
#' associations between loci (`genic`) and the remainder (`ld`), due to linkage
#' disequilibrium between loci and departures from Hardy-Weinberg. All
#' matrices use the `n - 1` denominator of [stats::cov()], so
#' `G = genic + ld`.
#'
#' * `method = "genic"` (exact; architectures where genotypic values are sums
#'   of allele effects, [additive()] and [binary_b()]): twice the sum over loci
#'   of the covariance of allele effects, pooling both alleles of every
#'   individual.
#' * `method = "shuffle"` (any architecture, including [multilinear()]): at
#'   each locus the `2n` allele copies are pooled and randomly reassigned to
#'   individuals, which removes associations between loci and within loci but
#'   keeps allele frequencies and every allele's full state (its effect vector,
#'   or its `y` and column of `B`). `genic` is the mean G of `shuffles` shuffled
#'   populations (do O & Whitlock 2023 call it the pleiotropic part of G). For
#'   additive architectures its expectation is the exact genic part.
#'
#' @param pop an `evomod_population`.
#' @param method `"genic"` or `"shuffle"`, see above.
#' @param shuffles number of shuffled populations averaged (`method = "shuffle"`).
#' @return A list with matrices `G`, `genic` and `ld`.
#' @references do O, I., & Whitlock, M. C. (2023). The evolution of genetic
#'   covariance and modularity as a result of multigenerational environmental
#'   fluctuation. *Evolution Letters*, 7(6), 457-466.
#' @examples
#' model <- evomod_model(additive(loci = 20, traits = 2),
#'                       mutation(rate = 0.01, cov = 0.1), init_sd = 0.1)
#' pop <- evolve(model, stabilizing(20, omega = 10), n = 100, seed = 1)
#' g_decomposition(pop)
#' g_decomposition(pop, method = "shuffle", shuffles = 20)
#' @export
g_decomposition <- function(pop, method = c("genic", "shuffle"), shuffles = 10) {
  method <- match.arg(method)
  arch <- pop$model$architecture
  G <- stats::cov(t(pop$x))
  if (method == "genic") {
    genic <- genic_covariance_cpp(genomes(pop), arch)
    if (!length(genic))
      stop("the genic part of G is only defined for architectures where genotypic values ",
           "are sums of allele effects (additive, binary_b); use method = \"shuffle\"",
           call. = FALSE)
  } else {
    check_count(shuffles, "shuffles")
    genic <- Reduce(`+`, lapply(seq_len(shuffles), function(i)
      stats::cov(t(genotypic_values(arch, shuffle_genomes(genomes(pop))))))) / shuffles
  }
  genic <- as.matrix(genic)
  nm <- trait_names(pop)
  dimnames(genic) <- dimnames(G) <- list(nm, nm)
  list(G = G, genic = genic, ld = G - genic)
}

#' Shuffle alleles among individuals
#'
#' At each locus, pools the `2n` allele copies of the population and assigns
#' them to individuals at random. Allele frequencies and the state of every
#' allele are kept; associations between loci (linkage disequilibrium) and
#' between the two alleles of a locus (departures from Hardy-Weinberg) are
#' removed.
#'
#' @param alleles genomes, an array `rows x alleles x individuals` (see
#'   [genomes()]).
#' @return An array of the same shape.
#' @examples
#' model <- evomod_model(additive(loci = 5, traits = 2), mutation(0.01, 0.1), init_sd = 0.1)
#' pop <- evolve(model, drift(5), n = 20, seed = 1)
#' shuffled <- shuffle_genomes(genomes(pop))
#' @export
shuffle_genomes <- function(alleles) {
  d <- dim(alleles)
  n <- d[3]
  out <- alleles
  for (l in seq_len(d[2] / 2)) {
    cols <- c(2 * l - 1, 2 * l)
    # columns: allele 1 of individual 1, allele 2 of individual 1, allele 1 of individual 2, ...
    pool <- matrix(alleles[, cols, , drop = FALSE], nrow = d[1])
    out[, cols, ] <- array(pool[, sample.int(2 * n), drop = FALSE], c(d[1], 2, n))
  }
  out
}

#' Mutational covariance
#'
#' Effects of new mutations on genotypic values in the current population.
#' Genomes are drawn at random (with replacement), and each receives exactly
#' one mutation event of the model's mutation process, at a site chosen in
#' proportion to the mutation rates; the effect is the change in its genotypic
#' value. With epistasis ([multilinear()]) or [binary_b()] the effects depend
#' on the genetic background, so the mutational covariance evolves.
#'
#' @param pop an `evomod_population`.
#' @param n number of mutations sampled.
#' @param effects if `TRUE`, also return the sampled effects.
#' @return A list with
#'   * `mean`: mean effect of a mutation on each trait (non-zero with
#'     directional epistasis);
#'   * `per_mutation`: covariance matrix of the effects of one mutation;
#'   * `expected_mutations`: expected number of mutation events per genome per
#'     generation;
#'   * `M`: mutational (co)variance input per generation,
#'     `expected_mutations * E[dx dx']`;
#'   * `effects`: the sampled effects (traits x n), if requested.
#' @examples
#' model <- evomod_model(additive(loci = 20, traits = 2),
#'                       mutation(rate = 0.01, cov = 0.1), init_sd = 0.1)
#' pop <- evolve(model, stabilizing(20, omega = 10), n = 100, seed = 1)
#' mutational_covariance(pop, n = 5000)
#' @export
mutational_covariance <- function(pop, n = 10000, effects = FALSE) {
  check_count(n, "n")
  model <- pop$model
  g <- genomes(pop)
  individuals <- sample.int(dim(g)[3], n, replace = TRUE) - 1L
  dx <- single_mutation_effects_cpp(g, model$traits, model$architecture, model$mutation,
                                    individuals)
  nm <- trait_names(pop)
  rownames(dx) <- nm
  second <- tcrossprod(dx) / n
  per_mutation <- stats::cov(t(dx))
  rate <- expected_mutations_cpp(matrix(g[, , 1], nrow = dim(g)[1]), model$mutation)
  dimnames(per_mutation) <- dimnames(second) <- list(nm, nm)
  out <- list(mean = rowMeans(dx), per_mutation = per_mutation,
              expected_mutations = rate, M = rate * second)
  if (effects) out$effects <- dx
  out
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
