#' Mutation models
#'
#' * `mutation()`: for [additive()] and [multilinear()] architectures. Each
#'   allele mutates with probability `rate` per generation; a mutation adds a
#'   draw from a multivariate normal with mean zero and covariance `cov` to the
#'   allele's effect vector.
#' * `binary_b_mutation()`: for [binary_b()]. Each allele's effect `y`
#'   mutates with probability `rate`, adding `Normal(0, sd)`; each entry of `B`
#'   mutates with probability `rate_b`. With `b_mutation = "random"` (the C
#'   code) a mutated entry gets a random value, 0 or 1, so it changes half of
#'   the time; with `"flip"` it always changes.
#'
#' @param rate mutation probability per allele per generation.
#' @param cov mutational covariance matrix (traits x traits), or a single
#'   number `s`, meaning `s^2` times the identity.
#' @param sd standard deviation of mutations in `y`.
#' @param rate_b mutation probability per entry of `B` per generation.
#' @param b_mutation `"random"` or `"flip"`, see above.
#' @return An `evomod_mutation` object.
#' @examples
#' mutation(rate = 0.005, cov = 0.05)
#' binary_b_mutation(rate = 0.0005, sd = 0.2, rate_b = 0.0001)
#' mutation(rate = 0.005, cov = matrix(c(1, 0.5, 0.5, 1), 2) * 0.01)
#' @export
mutation <- function(rate, cov) {
  if (!is.numeric(rate) || length(rate) != 1 || rate < 0 || rate > 1)
    stop("`rate` must be a probability", call. = FALSE)
  structure(list(type = "gaussian", rate = rate, cov = cov), class = "evomod_mutation")
}

#' @rdname mutation
#' @export
binary_b_mutation <- function(rate, sd, rate_b, b_mutation = c("random", "flip")) {
  b_mutation <- match.arg(b_mutation)
  for (r in list(rate, rate_b))
    if (!is.numeric(r) || length(r) != 1 || r < 0 || r > 1)
      stop("`rate` and `rate_b` must be probabilities", call. = FALSE)
  if (!is.numeric(sd) || length(sd) != 1 || sd < 0)
    stop("`sd` must be a non-negative number", call. = FALSE)
  structure(list(type = "binary_b", rate = rate, sd = sd, rate_b = rate_b,
                 b_mutation = b_mutation, random_state = b_mutation == "random"),
            class = "evomod_mutation")
}

#' Simulation model
#'
#' Combines a genetic architecture, a mutation model and environmental noise.
#'
#' @param architecture an architecture, see [architectures].
#' @param mutation a mutation model matching the architecture, see [mutation()].
#' @param env_sd standard deviation of environmental effects (`z = x + e`).
#' @param init_sd standard deviation of allelic effects in a new population
#'   (`y` for [binary_b()]); 0 (the default) starts every individual at the
#'   reference genotype.
#' @return An `evomod_model` object.
#' @examples
#' model <- evomod_model(additive(loci = 50, traits = 2),
#'                       mutation(rate = 0.005, cov = 0.05), env_sd = 1)
#' @export
evomod_model <- function(architecture, mutation, env_sd = 0, init_sd = 0) {
  if (!inherits(architecture, "evomod_architecture"))
    stop("`architecture` must be created with additive(), multilinear(), ...", call. = FALSE)
  if (!inherits(mutation, "evomod_mutation"))
    stop("`mutation` must be created with mutation()", call. = FALSE)
  p <- architecture$traits
  if (mutation$type != architecture$mutation_type)
    stop("the ", architecture$type, " architecture needs ",
         if (architecture$mutation_type == "binary_b") "binary_b_mutation()" else "mutation()",
         call. = FALSE)
  if (mutation$type == "gaussian") {
    cov <- mutation$cov
    if (length(cov) == 1) cov <- diag(cov^2, p)
    cov <- as.matrix(cov)
    if (!identical(dim(cov), c(p, p)) || !isSymmetric(unname(cov)))
      stop("mutational `cov` must be a symmetric ", p, " x ", p, " matrix", call. = FALSE)
    if (any(eigen(cov, symmetric = TRUE, only.values = TRUE)$values <= 0))
      stop("mutational `cov` must be positive definite", call. = FALSE)
    mutation$cov <- cov
  }
  structure(list(architecture = architecture, mutation = mutation,
                 env_sd = env_sd, init_sd = init_sd, traits = p,
                 loci = architecture$loci),
            class = "evomod_model")
}

#' @export
print.evomod_model <- function(x, ...) {
  cat("<evomod model>\n")
  cat("  architecture:", x$architecture$type, "-", x$loci, "loci,", x$traits, "traits\n")
  cat("  mutation:     rate", x$mutation$rate, "per allele\n")
  cat("  env_sd:      ", x$env_sd, "  init_sd:", x$init_sd, "\n")
  invisible(x)
}
