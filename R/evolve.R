#' Evolve a population
#'
#' Runs a population through a list of selection regimes. Starting from a
#' model creates a new population of `n` individuals; starting from a
#' population continues it (same individuals, generation count, optimum and
#' recorded statistics).
#'
#' Each generation: the optimum moves (directional regimes), fitness is
#' computed from the current phenotypes, every individual mutates, `n` pairs
#' of parents are drawn with probability proportional to fitness (with
#' replacement), and each offspring gets one random allele per locus from
#' each parent (free recombination). Random numbers come from R's generator,
#' so `set.seed()` (or `seed`) makes runs reproducible.
#'
#' @param x an `evomod_model` (new population) or an `evomod_population`
#'   (continue it).
#' @param regimes a regime or a list of regimes, see [regimes].
#' @param ... further arguments for the methods.
#' @return An `evomod_population`; see [trait_means()], [g_matrices()] and
#'   [variances()] for the recorded statistics.
#' @examples
#' model <- evomod_model(additive(loci = 20, traits = 2),
#'                       mutation(rate = 0.01, cov = 0.1), env_sd = 1, init_sd = 0.1)
#' pop <- evolve(model, n = 100, regimes = list(stabilizing(50, omega = 10)),
#'               record_every = 10, seed = 1)
#' pop <- evolve(pop, directional(50, delta_theta = c(0.02, 0), omega = 10))
#' trait_means(pop)
#' @export
evolve <- function(x, regimes, ...) UseMethod("evolve")

#' @rdname evolve
#' @param n population size.
#' @param record_every record population statistics every this many
#'   generations (and at the start).
#' @param seed optional seed passed to `set.seed()` before the run.
#' @param initial optional starting genomes, an array
#'   `rows x alleles x n` (see [genomes()]); by default they are drawn as
#'   described in [evomod_model()].
#' @export
evolve.evomod_model <- function(x, regimes, n, record_every = 1, seed = NULL,
                                initial = NULL, ...) {
  check_count(n, "n")
  if (!is.null(x$life_cycle) && x$life_cycle$type == "viability_selection" && n %% 2 != 0)
    stop("viability_selection() needs an even population size (half females)", call. = FALSE)
  check_count(record_every, "record_every")
  if (!is.null(seed)) set.seed(seed)
  p <- x$traits
  if (is.null(initial)) {
    alleles <- init_genomes(x$architecture, n, x$init_sd)
  } else {
    alleles <- initial
    if (!is.numeric(alleles) || length(dim(alleles)) != 3 ||
        !identical(as.integer(dim(alleles)), c(x$architecture$state_rows, 2L * x$loci, as.integer(n))))
      stop("`initial` must be a ", x$architecture$state_rows, " x ", 2 * x$loci, " x ", n,
           " array", call. = FALSE)
    storage.mode(alleles) <- "double"
  }
  pop <- structure(list(model = x, n = as.integer(n), alleles = alleles,
                        x = NULL, z = NULL, generation = 0L, theta = rep(0, p),
                        record_every = as.integer(record_every),
                        records = NULL, regimes = list()),
                   class = "evomod_population")
  run_regimes(pop, regimes, record_start = TRUE)
}

#' @rdname evolve
#' @export
evolve.evomod_population <- function(x, regimes, record_every = x$record_every,
                                     seed = NULL, ...) {
  check_count(record_every, "record_every")
  if (!is.null(seed)) set.seed(seed)
  x$record_every <- as.integer(record_every)
  run_regimes(x, regimes, record_start = FALSE)
}

run_regimes <- function(pop, regimes, record_start) {
  model <- pop$model
  resolved <- resolve_regimes(regimes, model$traits, pop$theta)
  z <- if (is.null(pop$z)) matrix(0, 0, 0) else pop$z
  out <- evolve_cpp(pop$alleles, model$traits, z, model$architecture, model$mutation,
                    model$recombination, model$env_sd,
                    resolved, life_cycle_of(model), pop$generation, pop$record_every,
                    record_start)
  pop$alleles <- out$alleles
  pop$x <- out$x
  pop$z <- out$z
  pop$theta <- as.numeric(out$theta)
  pop$generation <- out$generation
  pop$records <- append_records(pop$records, out$records)
  start <- pop$generation - sum(vapply(resolved, `[[`, 1L, "generations"))
  pop$regimes <- c(pop$regimes, lapply(resolved, function(r) {
    r$start <- start; start <<- start + r$generations; r
  }))
  pop
}

life_cycle_of <- function(model) if (is.null(model$life_cycle)) wright_fisher() else model$life_cycle

append_records <- function(old, new) {
  new$generation <- as.integer(new$generation)
  if (is.null(old)) return(new)
  bind3 <- function(a, b) {
    out <- array(0, c(dim(a)[1:2], dim(a)[3] + dim(b)[3]))
    out[, , seq_len(dim(a)[3])] <- a
    out[, , dim(a)[3] + seq_len(dim(b)[3])] <- b
    out
  }
  out <- list(generation = c(old$generation, new$generation),
              theta = cbind(old$theta, new$theta),
              mean_x = cbind(old$mean_x, new$mean_x),
              mean_z = cbind(old$mean_z, new$mean_z),
              G = bind3(old$G, new$G), P = bind3(old$P, new$P),
       mean_fitness = c(old$mean_fitness, new$mean_fitness))
  if (!is.null(old$G_genic) && !is.null(new$G_genic)) out$G_genic <- bind3(old$G_genic, new$G_genic)
  out
}

#' @export
print.evomod_population <- function(x, ...) {
  cat("<evomod population>", x$n, "individuals at generation", x$generation, "\n")
  cat("  architecture:", x$model$architecture$type, "-", x$model$loci, "loci,",
      x$model$traits, "traits\n")
  cat("  regimes:", paste0(vapply(x$regimes, `[[`, "", "type"), " (",
                           vapply(x$regimes, `[[`, 1L, "generations"), ")",
                           collapse = ", "), "\n")
  cat("  records:", length(x$records$generation), "generations\n")
  invisible(x)
}
