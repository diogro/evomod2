#' Genetic architectures
#'
#' An architecture maps a genome to genotypic values. Genomes are matrices of
#' allelic effect vectors (traits x alleles); alleles \eqn{2l - 1}{2l - 1} and
#' \eqn{2l}{2l} belong to locus \eqn{l}{l}, and the effect of a locus is the sum
#' of its two alleles (\eqn{Y_l}{Y_l}). Effects are measured from the reference genotype, where all
#' effects are zero.
#'
#' * `additive()`: \eqn{x = \sum_l Y_l}{x = sum_l Y_l}.
#' * `binary_b()`: the binary ontogenetic matrix model of Melo & Marroig
#'   (2015). Each allele has a scalar effect \eqn{y}{y} and a column of a binary
#'   matrix \eqn{B}{B} saying which traits it affects; \eqn{x = B y}{x = B y}. Use with
#'   [binary_b_mutation()]. Genomes have `traits + 1` rows: \eqn{y}{y}, then \eqn{B}{B}.
#' * `contribution_matrix()`: the model of do O & Whitlock (2023). `genes`
#'   weight genes \eqn{a_j}{a_j} and a `traits x genes` matrix of contribution
#'   genes \eqn{c_{ij}}{c_ij}, each its own locus with a scalar value (the two alleles of a
#'   locus are summed); \eqn{x_i = \sum_j c_{ij} a_j}{x_i = sum_j c_ij a_j}. Pleiotropy
#'   (\eqn{C}{C}) evolves
#'   continuously and can be negative, and the map is bilinear (epistasis
#'   between weight and contribution genes). Loci `1..genes` are the weight
#'   genes; locus `genes + (j - 1) * traits + i` is \eqn{c_{ij}}{c_ij}. Use with
#'   [contribution_mutation()]. New populations start with every individual
#'   at `init_weights` and `init_contributions` (each allele half of the
#'   locus value), plus `Normal(0, init_sd)` noise on every allele.
#' * `multilinear()` with a pairwise `epsilon` ([pairwise_epsilon()]): a
#'   different tensor for every pair of loci (Jones, Burger & Arnold 2014),
#'   \eqn{x = S + \sum_{i<j} \varepsilon_{ij}(Y_i, Y_j)}{x = S + sum_{i<j} eps_ij(Y_i, Y_j)}.
#' * `multilinear()`: pairwise multilinear epistasis (Hansen & Wagner 2001),
#'   \eqn{x = S + \sum_{i<j} \varepsilon(Y_i, Y_j)}{x = S + sum_{i<j} eps(Y_i, Y_j)},
#'   with \eqn{S = \sum_l Y_l}{S = sum_l Y_l} and
#'   \eqn{\varepsilon(a, b)_k = \sum_{l,m} \varepsilon_{klm} a_l b_m}{eps(a, b)[k] = sum_{l,m} epsilon[k, l, m] a[l] b[m]}. The same tensor is
#'   used for every pair of loci and must be symmetric in its last two indices.
#'
#' @param loci number of diploid loci.
#' @param traits number of traits.
#' @param epsilon epistasis coefficients: a `traits x traits x traits` array
#'   shared by all pairs of loci (symmetric in its last two indices; see
#'   [random_epsilon()]), or a `traits x traits x traits x loci x loci` array
#'   with a different tensor for every pair of loci, as in Jones, Burger &
#'   Arnold (2014) (see [pairwise_epsilon()]). For the pairwise form,
#'   `epsilon[k, a, b, i, j]` multiplies `Y_i[a] Y_j[b]` in trait `k`, and only
#'   entries with `i < j` are used.
#' @param genes number of weight genes.
#' @param init_weights starting values of the weight genes (one per gene).
#' @param init_contributions starting contribution matrix (`traits x genes`);
#'   by default ones on the diagonal (gene `j` affects trait `j` only, the
#'   identity when `genes == traits`).
#' @param init_b probability that an entry of `B` is 1 in a new population
#'   (0.25 in the C code, which starts from random mutations of an empty `B`).
#' @param init_y how allelic effects `y` of a new population are drawn, with
#'   the model's `init_sd`: `"normal"`, `Normal(0, init_sd)` for every allele;
#'   `"mutations"`, as in the C code, `2 x loci` mutations of
#'   `Normal(0, init_sd)` added to alleles chosen at random in each genome.
#' @return An `evomod_architecture` object.
#' @examples
#' additive(loci = 50, traits = 2)
#' binary_b(loci = 50, traits = 4)
#' contribution_matrix(traits = 4)
#' multilinear(loci = 50, traits = 2, epsilon = random_epsilon(2, sd = 0.05))
#' @name architectures
NULL

#' @rdname architectures
#' @export
additive <- function(loci, traits) {
  check_count(loci, "loci")
  check_count(traits, "traits")
  structure(list(type = "additive", loci = as.integer(loci), traits = as.integer(traits),
                 state_rows = as.integer(traits), mutation_type = "gaussian"),
            class = "evomod_architecture")
}

#' @rdname architectures
#' @export
binary_b <- function(loci, traits, init_b = 0.25, init_y = c("normal", "mutations")) {
  init_y <- match.arg(init_y)
  check_count(loci, "loci")
  check_count(traits, "traits")
  if (!is.numeric(init_b) || length(init_b) != 1 || init_b < 0 || init_b > 1)
    stop("`init_b` must be a probability", call. = FALSE)
  structure(list(type = "binary_b", loci = as.integer(loci), traits = as.integer(traits),
                 state_rows = as.integer(traits) + 1L, mutation_type = "binary_b",
                 init_b = init_b, init_y = init_y),
            class = "evomod_architecture")
}

#' @rdname architectures
#' @export
contribution_matrix <- function(traits, genes = traits, init_weights = rep(1, genes),
                                init_contributions = diag(1, traits, genes)) {
  check_count(traits, "traits")
  check_count(genes, "genes")
  p <- as.integer(traits); k <- as.integer(genes)
  if (!is.numeric(init_weights) || length(init_weights) != k)
    stop("`init_weights` must have one value per gene (", k, ")", call. = FALSE)
  init_contributions <- as.matrix(init_contributions)
  if (!identical(dim(init_contributions), c(p, k)))
    stop("`init_contributions` must be a ", p, " x ", k, " matrix", call. = FALSE)
  structure(list(type = "contribution_matrix", loci = k + p * k, traits = p, genes = k,
                 state_rows = 1L, mutation_type = "contribution",
                 init_weights = as.numeric(init_weights),
                 init_contributions = init_contributions),
            class = "evomod_architecture")
}

#' @rdname architectures
#' @export
multilinear <- function(loci, traits, epsilon) {
  check_count(loci, "loci")
  check_count(traits, "traits")
  p <- as.integer(traits)
  L <- as.integer(loci)
  epsilon <- as.array(epsilon)
  if (length(dim(epsilon)) == 5) {
    if (!is.numeric(epsilon) || !identical(dim(epsilon), c(p, p, p, L, L)))
      stop("pairwise `epsilon` must be a numeric ", p, " x ", p, " x ", p, " x ", L, " x ", L,
           " array", call. = FALSE)
    pairs <- which(upper.tri(diag(L)), arr.ind = TRUE)
    pairs <- pairs[order(pairs[, 1], pairs[, 2]), , drop = FALSE]   # (1,2), (1,3), ..., (2,3), ...
    slices <- array(0, c(p, p, p * nrow(pairs)))
    for (q in seq_len(nrow(pairs))) for (k in seq_len(p))
      slices[, , (q - 1) * p + k] <- epsilon[k, , , pairs[q, 1], pairs[q, 2]]
    return(structure(list(type = "multilinear_pairwise", loci = L, traits = p, state_rows = p,
                          mutation_type = "gaussian", epsilon = epsilon, pair_slices = slices),
                     class = "evomod_architecture"))
  }
  if (!is.numeric(epsilon) || !identical(dim(epsilon), c(p, p, p)))
    stop("`epsilon` must be a numeric ", p, " x ", p, " x ", p, " array")
  for (k in seq_len(p))
    if (!isSymmetric(unname(matrix(epsilon[k, , ], p, p)), tol = 1e-10))
      stop("`epsilon[", k, ", , ]` must be symmetric")
  structure(list(type = "multilinear", loci = as.integer(loci), traits = p,
                 state_rows = p, mutation_type = "gaussian", epsilon = epsilon,
                 # slice k (third index) is epsilon[k, , ], as the C++ code expects
                 epsilon_slices = aperm(epsilon, c(2, 3, 1))),
            class = "evomod_architecture")
}

#' Random pairwise epistasis coefficients
#'
#' A different epistasis tensor for every pair of loci, with entries drawn
#' from `Normal(mean, sd)`, as in Jones, Burger & Arnold (2014), where the
#' variance of the coefficients is `sd^2`. With `within_trait = TRUE` only the
#' coefficients for interactions of a trait with itself (`epsilon[k, k, k, i, j]`)
#' are non-zero; `traits_with_epistasis` restricts epistasis to some traits.
#'
#' @param loci number of loci.
#' @param traits number of traits.
#' @param sd,mean standard deviation and mean of the coefficients.
#' @param within_trait only interactions of each trait with itself.
#' @param traits_with_epistasis traits whose values have epistatic terms.
#' @return A `traits x traits x traits x loci x loci` array (zero for `i >= j`),
#'   for [multilinear()].
#' @examples
#' eps <- pairwise_epsilon(loci = 10, traits = 2, sd = 1)
#' multilinear(loci = 10, traits = 2, epsilon = eps)
#' @export
pairwise_epsilon <- function(loci, traits, sd = 1, mean = 0, within_trait = FALSE,
                             traits_with_epistasis = seq_len(traits)) {
  p <- as.integer(traits); L <- as.integer(loci)
  e <- array(0, c(p, p, p, L, L))
  for (i in seq_len(L - 1)) for (j in (i + 1):L) {
    t <- array(stats::rnorm(p^3, mean, sd), c(p, p, p))
    if (within_trait) {
      keep <- array(0, c(p, p, p))
      for (k in seq_len(p)) keep[k, k, k] <- 1
      t <- t * keep
    }
    t[setdiff(seq_len(p), traits_with_epistasis), , ] <- 0
    e[, , , i, j] <- t
  }
  e
}

#' Random epistasis tensor
#'
#' Entries are drawn from `Normal(mean, sd)` and symmetrised in the last two
#' indices. A positive `mean` gives positive directional epistasis on average.
#'
#' @param traits number of traits.
#' @param sd,mean standard deviation and mean of the entries.
#' @return A `traits x traits x traits` array.
#' @examples
#' random_epsilon(2, sd = 0.1)
#' @export
random_epsilon <- function(traits, sd = 1, mean = 0) {
  p <- as.integer(traits)
  e <- array(stats::rnorm(p^3, mean, sd), c(p, p, p))
  for (k in seq_len(p)) e[k, , ] <- (e[k, , ] + t(e[k, , ])) / 2
  e
}

#' @export
print.evomod_architecture <- function(x, ...) {
  cat("<evomod architecture>", x$type, "-", x$loci, "loci,", x$traits, "traits\n")
  invisible(x)
}

# Genomes of a new population: genome rows x alleles x individuals.
init_genomes <- function(arch, n, init_sd) {
  n_alleles <- 2L * arch$loci
  if (arch$type == "contribution_matrix") {
    values <- c(arch$init_weights, as.vector(arch$init_contributions))   # c_ij column-major
    g <- array(rep(rep(values / 2, each = 2), n), c(1, n_alleles, n))
    return(g + stats::rnorm(n_alleles * n, 0, init_sd))
  }
  if (arch$type == "binary_b") {
    g <- array(0, c(arch$state_rows, n_alleles, n))
    if (arch$init_y == "normal") {
      g[1, , ] <- stats::rnorm(n_alleles * n, 0, init_sd)
    } else {
      for (i in seq_len(n)) {
        hits <- tabulate(sample.int(n_alleles, n_alleles, replace = TRUE), n_alleles)
        g[1, , i] <- stats::rnorm(n_alleles, 0, init_sd * sqrt(hits))
      }
    }
    g[-1, , ] <- stats::rbinom(arch$traits * n_alleles * n, 1, arch$init_b)
    return(g)
  }
  array(stats::rnorm(arch$state_rows * n_alleles * n, 0, init_sd),
        c(arch$state_rows, n_alleles, n))
}

check_count <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < 1 || x != round(x))
    stop("`", name, "` must be a positive whole number", call. = FALSE)
}
