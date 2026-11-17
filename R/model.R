#' Mutation models
#'
#' * `mutation()`: for [additive()] and [multilinear()] architectures;
#'   `hoc_mutation()` is the house of cards version (`persistence = 0`). Each
#'   allele mutates with probability `rate` per generation; a mutation adds a
#'   draw from a multivariate normal with mean zero and covariance `cov` to the
#'   allele's effect vector.
#' * `contribution_mutation()`: for [contribution_matrix()]. Each allele of
#'   a weight gene mutates with probability `rate_weight`, each allele of a
#'   contribution gene with `rate_contribution`; a mutation adds
#'   `Normal(0, sd)` to the allele's value.
#' * `binary_b_mutation()`: for [binary_b()]. Each allele's effect `y`
#'   mutates with probability `rate`, adding `Normal(0, sd)`; each entry of `B`
#'   mutates with probability `rate_b`. With `b_mutation = "flip"` (the
#'   default, and the current C code) a mutated entry always changes; with
#'   `"random"` it gets a random value, 0 or 1, so it changes half of the time,
#'   as in the C code used for the simulations of Melo & Marroig (2015).
#'
#' @param rate mutation probability per allele per generation.
#' @param cov mutational covariance matrix (traits x traits), or a single
#'   number \eqn{s}{s}, meaning \eqn{s^2 I}{s^2 times the identity}.
#' @param persistence how much of the old effect a mutation keeps: the new
#'   effect is \eqn{\rho e + \delta}{persistence * old + draw}, with \eqn{\rho}{rho} =
#'   `persistence`, \eqn{e}{e} the old effect and \eqn{\delta}{delta} a fresh draw. 1 (the default) is a random walk,
#'   where mutational effects accumulate; 0 is the house of cards model
#'   (Kingman 1978; Turelli 1984), where a mutation replaces the effect with a
#'   fresh draw unrelated to it; values in between give a regression towards
#'   the reference genotype (Zeng & Cockerham 1993). For [binary_b()] it applies
#'   to `y`.
#' @param sd standard deviation of mutations in `y`.
#' @param rate_b mutation probability per entry of `B` per generation.
#' @param b_mutation `"random"` or `"flip"`, see above.
#' @return An `evomod_mutation` object.
#' @examples
#' mutation(rate = 0.005, cov = 0.05)
#' hoc_mutation(rate = 0.005, cov = 0.5)   # house of cards
#' binary_b_mutation(rate = 0.0005, sd = 0.2, rate_b = 0.0001)
#' mutation(rate = 0.005, cov = matrix(c(1, 0.5, 0.5, 1), 2) * 0.01)
#' @export
mutation <- function(rate, cov, persistence = 1) {
  if (!is.numeric(rate) || length(rate) != 1 || rate < 0 || rate > 1)
    stop("`rate` must be a probability", call. = FALSE)
  check_persistence(persistence)
  structure(list(type = "gaussian", rate = rate, cov = cov, persistence = persistence),
            class = "evomod_mutation")
}

#' @rdname mutation
#' @export
hoc_mutation <- function(rate, cov) mutation(rate, cov, persistence = 0)

check_persistence <- function(persistence) {
  if (!is.numeric(persistence) || length(persistence) != 1 || persistence < 0 || persistence > 1)
    stop("`persistence` must be a number between 0 and 1", call. = FALSE)
}

#' @rdname mutation
#' @param rate_weight,rate_contribution mutation probability per allele per
#'   generation of weight genes and of contribution genes.
#' @export
contribution_mutation <- function(rate_weight, rate_contribution, sd, persistence = 1) {
  check_persistence(persistence)
  for (r in list(rate_weight, rate_contribution))
    if (!is.numeric(r) || length(r) != 1 || r < 0 || r > 1)
      stop("mutation rates must be probabilities", call. = FALSE)
  if (!is.numeric(sd) || length(sd) != 1 || sd < 0)
    stop("`sd` must be a non-negative number", call. = FALSE)
  structure(list(type = "contribution", rate_weight = rate_weight,
                 rate_contribution = rate_contribution, sd = sd, persistence = persistence),
            class = "evomod_mutation")
}

#' @rdname mutation
#' @export
binary_b_mutation <- function(rate, sd, rate_b, b_mutation = c("flip", "random"),
                              persistence = 1) {
  b_mutation <- match.arg(b_mutation)
  check_persistence(persistence)
  for (r in list(rate, rate_b))
    if (!is.numeric(r) || length(r) != 1 || r < 0 || r > 1)
      stop("`rate` and `rate_b` must be probabilities", call. = FALSE)
  if (!is.numeric(sd) || length(sd) != 1 || sd < 0)
    stop("`sd` must be a non-negative number", call. = FALSE)
  structure(list(type = "binary_b", rate = rate, sd = sd, rate_b = rate_b,
                 b_mutation = b_mutation, random_state = b_mutation == "random",
                 persistence = persistence),
            class = "evomod_mutation")
}

#' Life cycles
#'
#' * `wright_fisher()` (the default): discrete generations. Every generation,
#'   each individual mutates and the whole population is replaced by `n`
#'   offspring of parents drawn with probability proportional to fitness.
#' * `moran()`: overlapping generations. A generation is `n` birth-death
#'   events; at each event a random individual dies and is replaced by the
#'   offspring of two parents drawn with probability proportional to fitness
#'   (they can be the same individual, as in `wright_fisher()`). The offspring
#'   mutates with the model's per-generation rates, so the mutational input per
#'   generation is the same as with `wright_fisher()`. Drift per generation is
#'   about twice as strong as with `wright_fisher()` at the same `n`.
#'
#' * `viability_selection()`: the life cycle of Jones, Arnold & Burger
#'   (2003, 2004, 2014). Half the population is female; each female mates with
#'   `mates` males drawn at random (males can mate many times) and has
#'   `offspring_per_mate` offspring with each. Zygotes mutate, survive with
#'   probability equal to their fitness (absolute, `W(z) <= 1`), and `n`
#'   adults are drawn at random from the survivors (if fewer than `n` survive,
#'   survivors are resampled to keep `n`). `n` must be even. With the defaults
#'   (2 mates, 2 offspring each), the effective population size is about `n`,
#'   as with `wright_fisher()`.
#'
#' Statistics are recorded, and directional or fluctuating optima move,
#' between generations.
#'
#' @param mates number of males each female mates with.
#' @param offspring_per_mate number of offspring per female and male.
#'
#' @return An `evomod_life_cycle` object, for [evomod_model()].
#' @examples
#' moran()
#' @name life_cycles
NULL

#' @rdname life_cycles
#' @export
wright_fisher <- function() structure(list(type = "wright_fisher"), class = "evomod_life_cycle")

#' @rdname life_cycles
#' @export
moran <- function() structure(list(type = "moran"), class = "evomod_life_cycle")

#' @rdname life_cycles
#' @export
viability_selection <- function(mates = 2, offspring_per_mate = 2) {
  check_count(mates, "mates")
  check_count(offspring_per_mate, "offspring_per_mate")
  if (mates * offspring_per_mate < 2)
    stop("each female needs at least 2 offspring (mates x offspring_per_mate)", call. = FALSE)
  structure(list(type = "viability_selection", mates = as.integer(mates),
                 offspring_per_mate = as.integer(offspring_per_mate)),
            class = "evomod_life_cycle")
}

#' Simulation model
#'
#' Combines a genetic architecture, a mutation model, recombination and
#' environmental noise.
#'
#' @param architecture an architecture, see [architectures].
#' @param mutation a mutation model matching the architecture, see [mutation()].
#' @param recombination how loci are inherited, see [recombination]; free
#'   recombination by default.
#' @param life_cycle [wright_fisher()] (discrete generations, the default),
#'   [moran()] (overlapping generations) or [viability_selection()] (separate
#'   sexes and viability selection, Jones et al.).
#' @param env_sd standard deviation of environmental effects \eqn{e}{e} in
#'   \eqn{z = x + e}{z = x + e}.
#' @param init_sd standard deviation of allelic effects in a new population
#'   (`y` for [binary_b()]); 0 (the default) starts every individual at the
#'   reference genotype.
#' @return An `evomod_model` object.
#' @examples
#' model <- evomod_model(additive(loci = 50, traits = 2),
#'                       mutation(rate = 0.005, cov = 0.05), env_sd = 1)
#' @export
evomod_model <- function(architecture, mutation, env_sd = 0, init_sd = 0,
                         recombination = free_recombination(), life_cycle = wright_fisher()) {
  if (!inherits(life_cycle, "evomod_life_cycle"))
    stop("`life_cycle` must be wright_fisher() or moran()", call. = FALSE)
  if (!inherits(architecture, "evomod_architecture"))
    stop("`architecture` must be created with additive(), multilinear(), ...", call. = FALSE)
  if (!inherits(mutation, "evomod_mutation"))
    stop("`mutation` must be created with mutation()", call. = FALSE)
  p <- architecture$traits
  if (mutation$type != architecture$mutation_type)
    stop("the ", architecture$type, " architecture needs ",
         switch(architecture$mutation_type, binary_b = "binary_b_mutation()",
                contribution = "contribution_mutation()", "mutation()"),
         call. = FALSE)
  if (mutation$type == "contribution") mutation$weight_loci <- architecture$genes
  if (mutation$type == "gaussian") {
    cov <- mutation$cov
    if (length(cov) == 1 && !is.matrix(cov)) cov <- diag(cov^2, p)
    cov <- as.matrix(cov)
    if (!identical(dim(cov), c(p, p)) || !isSymmetric(unname(cov)))
      stop("mutational `cov` must be a symmetric ", p, " x ", p, " matrix", call. = FALSE)
    if (any(eigen(cov, symmetric = TRUE, only.values = TRUE)$values <= 0))
      stop("mutational `cov` must be positive definite", call. = FALSE)
    mutation$cov <- cov
  }
  recombination <- resolve_recombination(recombination, architecture$loci)
  structure(list(architecture = architecture, mutation = mutation,
                 recombination = recombination, life_cycle = life_cycle,
                 env_sd = env_sd, init_sd = init_sd, traits = p,
                 loci = architecture$loci),
            class = "evomod_model")
}

#' @export
print.evomod_model <- function(x, ...) {
  cat("<evomod model>\n")
  cat("  architecture:", x$architecture$type, "-", x$loci, "loci,", x$traits, "traits\n")
  pers <- if (is.null(x$mutation$persistence)) 1 else x$mutation$persistence
  cat("  mutation:     rate", x$mutation$rate, "per allele,",
      if (pers == 1) "random walk" else if (pers == 0) "house of cards" else
        paste("persistence", pers), "\n")
  r <- x$recombination
  cat("  recombination:", if (r$type == "free") "free" else
      if (r$map_type == "genetic_map") paste(r$chromosomes, "chromosome(s) of", r$length, "Morgans")
      else "linkage map", "\n")
  cat("  life cycle:  ", if (is.null(x$life_cycle)) "wright_fisher" else x$life_cycle$type, "\n")
  cat("  env_sd:      ", x$env_sd, "  init_sd:", x$init_sd, "\n")
  invisible(x)
}
