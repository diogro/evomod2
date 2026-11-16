#' Genetic architectures
#'
#' An architecture maps a genome to genotypic values. Genomes are matrices of
#' allelic effect vectors (traits x alleles); alleles `2l - 1` and `2l` belong
#' to locus `l`, and the effect of a locus is the sum of its two alleles
#' (`Y_l`). Effects are measured from the reference genotype, where all
#' effects are zero.
#'
#' * `additive()`: `x = sum_l Y_l`.
#' * `binary_b()`: the binary ontogenetic matrix model of Melo & Marroig
#'   (2015). Each allele has a scalar effect `y` and a column of a binary
#'   matrix `B` saying which traits it affects; `x = B y`. Use with
#'   [binary_b_mutation()]. Genomes have `traits + 1` rows: `y`, then `B`.
#' * `multilinear()`: pairwise multilinear epistasis (Hansen & Wagner 2001),
#'   `x = S + sum_{i<j} eps(Y_i, Y_j)`, with `S = sum_l Y_l` and
#'   `eps(a, b)[k] = sum_{l,m} epsilon[k, l, m] a[l] b[m]`. The same tensor is
#'   used for every pair of loci and must be symmetric in its last two indices.
#'
#' @param loci number of diploid loci.
#' @param traits number of traits.
#' @param epsilon epistasis tensor, a `traits x traits x traits` array
#'   (see [random_epsilon()]).
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
multilinear <- function(loci, traits, epsilon) {
  check_count(loci, "loci")
  check_count(traits, "traits")
  p <- as.integer(traits)
  epsilon <- as.array(epsilon)
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
