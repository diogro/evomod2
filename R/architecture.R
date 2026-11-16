#' Genetic architectures
#'
#' An architecture maps a genome to genotypic values. Genomes are matrices of
#' allelic effect vectors (traits x alleles); alleles `2l - 1` and `2l` belong
#' to locus `l`, and the effect of a locus is the sum of its two alleles
#' (`Y_l`). Effects are measured from the reference genotype, where all
#' effects are zero.
#'
#' * `additive()`: `x = sum_l Y_l`.
#' * `multilinear()`: pairwise multilinear epistasis (Hansen & Wagner 2001),
#'   `x = S + sum_{i<j} eps(Y_i, Y_j)`, with `S = sum_l Y_l` and
#'   `eps(a, b)[k] = sum_{l,m} epsilon[k, l, m] a[l] b[m]`. The same tensor is
#'   used for every pair of loci and must be symmetric in its last two indices.
#'
#' @param loci number of diploid loci.
#' @param traits number of traits.
#' @param epsilon epistasis tensor, a `traits x traits x traits` array
#'   (see [random_epsilon()]).
#' @return An `evomod_architecture` object.
#' @examples
#' additive(loci = 50, traits = 2)
#' multilinear(loci = 50, traits = 2, epsilon = random_epsilon(2, sd = 0.05))
#' @name architectures
NULL

#' @rdname architectures
#' @export
additive <- function(loci, traits) {
  check_count(loci, "loci")
  check_count(traits, "traits")
  structure(list(type = "additive", loci = as.integer(loci), traits = as.integer(traits)),
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
                 epsilon = epsilon,
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

check_count <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < 1 || x != round(x))
    stop("`", name, "` must be a positive whole number", call. = FALSE)
}
