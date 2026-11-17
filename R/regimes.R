#' Selection regimes
#'
#' A simulation runs through a list of regimes:
#'
#' * `drift()`: no selection.
#' * `stabilizing()`: Gaussian fitness surface with optimum `theta` and
#'   covariance `omega`,
#'   \deqn{w(z) = \exp\left(-\frac{1}{2} (z - \theta)^\top \Omega^{-1} (z - \theta)\right).}{w(z) = exp(-(z - theta)' omega^-1 (z - theta) / 2).}
#' * `directional()`: as stabilizing, but the optimum moves by `delta_theta`
#'   every generation (before selection).
#' * `fluctuating()`: as stabilizing, but after every `every` generations the
#'   optimum jumps by `step` up or down (do O & Whitlock 2023). With
#'   `structure = "correlated"` all traits jump in the same direction; with
#'   `"independent"` each trait draws its own direction; a list of trait
#'   blocks (e.g. `list(1:2, 3:4)`) makes traits jump together within blocks
#'   and independently between them (traits not in any block jump alone).
#'   A jump that would take a trait of the block outside `bounds` is not
#'   taken, so at a bound the optimum stays or moves inward with probability
#'   1/2 each; `bounds = NULL` leaves the optimum unbounded.
#'
#' When `theta` is `NULL` the regime starts from the optimum at the end of the
#' previous regime (the origin for the first one).
#'
#' @param generations number of generations.
#' @param omega fitness surface covariance matrix (traits x traits), or a
#'   single number `w`, meaning `w` times the identity.
#' @param theta optimum at the start of the regime, or `NULL`.
#' @param delta_theta change in the optimum per generation (vector, one value
#'   per trait).
#' @param every number of generations between jumps of the optimum.
#' @param step size of a jump (one value, or one per trait).
#' @param structure `"correlated"`, `"independent"` or a list of trait blocks.
#' @param bounds lower and upper bound of the optimum (`NULL` for none).
#' @return An `evomod_regime` object.
#' @examples
#' regimes <- list(drift(100), stabilizing(1000, omega = 20),
#'                 directional(1000, delta_theta = c(0.005, 0), omega = 20))
#' fluctuating(5000, every = 10, step = 0.2, omega = 1, structure = "correlated")
#' fluctuating(5000, every = 10, step = 0.2, omega = 1, structure = list(1:2, 3:4))
#' @name regimes
NULL

new_regime <- function(type, generations, omega = NULL, theta = NULL, delta_theta = NULL, ...) {
  check_count(generations, "generations")
  structure(list(type = type, generations = as.integer(generations), omega = omega,
                 theta = theta, delta_theta = delta_theta, ...),
            class = "evomod_regime")
}

#' @rdname regimes
#' @export
drift <- function(generations) new_regime("drift", generations)

#' @rdname regimes
#' @export
stabilizing <- function(generations, omega, theta = NULL)
  new_regime("stabilizing", generations, omega, theta)

#' @rdname regimes
#' @export
directional <- function(generations, delta_theta, omega, theta = NULL)
  new_regime("directional", generations, omega, theta, delta_theta)

#' @rdname regimes
#' @export
fluctuating <- function(generations, every, step, omega, structure = "correlated",
                        bounds = c(-1, 1), theta = NULL) {
  check_count(every, "every")
  if (!is.numeric(step) || any(step < 0)) stop("`step` must be non-negative", call. = FALSE)
  if (!is.null(bounds) && (!is.numeric(bounds) || length(bounds) != 2 || bounds[1] >= bounds[2]))
    stop("`bounds` must be c(lower, upper) with lower < upper, or NULL", call. = FALSE)
  if (!(is.list(structure) || structure %in% c("correlated", "independent")))
    stop("`structure` must be \"correlated\", \"independent\" or a list of trait blocks",
         call. = FALSE)
  new_regime("fluctuating", generations, omega, theta, every = as.integer(every), step = step,
             structure = structure, bounds = bounds)
}

#' @export
print.evomod_regime <- function(x, ...) {
  cat("<evomod regime>", x$type, "for", x$generations, "generations\n")
  invisible(x)
}

# Regimes as the C++ code expects them, with the optimum at the start of each
# regime chained from the previous one.
resolve_regimes <- function(regimes, p, theta) {
  first <- TRUE
  if (inherits(regimes, "evomod_regime")) regimes <- list(regimes)
  if (!length(regimes) || !all(vapply(regimes, inherits, TRUE, "evomod_regime")))
    stop("`regimes` must be a list of drift(), stabilizing(), directional() or fluctuating() regimes",
         call. = FALSE)
  lapply(regimes, function(r) {
    omega <- if (is.null(r$omega)) diag(p) else r$omega
    if (length(omega) == 1 && !is.matrix(omega)) omega <- diag(omega, p)
    omega <- as.matrix(omega)
    if (!identical(dim(omega), c(p, p)) || !isSymmetric(unname(omega)))
      stop("`omega` must be a symmetric ", p, " x ", p, " matrix", call. = FALSE)
    reset <- first || !is.null(r$theta)
    if (!is.null(r$theta)) theta <- check_vector(r$theta, p, "theta")
    first <<- FALSE
    delta <- if (is.null(r$delta_theta)) rep(0, p) else check_vector(r$delta_theta, p, "delta_theta")
    out <- list(generations = r$generations, selective = r$type != "drift",
                reset_theta = reset, theta = theta, delta_theta = delta, omega = omega,
                type = r$type, every = 0L)
    if (r$type == "fluctuating") {
      step <- if (length(r$step) == 1) rep(r$step, p) else check_vector(r$step, p, "step")
      block <- fluctuation_blocks(r$structure, p)
      bounds <- if (is.null(r$bounds)) c(-Inf, Inf) else r$bounds
      if (any(theta < bounds[1] | theta > bounds[2]))
        stop("the optimum at the start of a fluctuating regime is outside `bounds`", call. = FALSE)
      out[c("every", "step", "block", "lower", "upper")] <-
        list(r$every, step, block, rep(bounds[1], p), rep(bounds[2], p))
    }
    out
  })
}

# 0-based block of each trait
fluctuation_blocks <- function(structure, p) {
  if (identical(structure, "correlated")) return(rep(0L, p))
  if (identical(structure, "independent")) return(seq_len(p) - 1L)
  traits <- unlist(structure)
  if (!is.numeric(traits) || any(traits < 1 | traits > p | traits != round(traits)) ||
      anyDuplicated(traits))
    stop("blocks must be disjoint sets of trait indices between 1 and ", p, call. = FALSE)
  block <- integer(p)
  for (b in seq_along(structure)) block[structure[[b]]] <- b - 1L
  alone <- setdiff(seq_len(p), traits)
  block[alone] <- length(structure) + seq_along(alone) - 1L
  block
}

check_vector <- function(x, p, name) {
  if (!is.numeric(x) || length(x) != p)
    stop("`", name, "` must be a numeric vector of length ", p, call. = FALSE)
  as.numeric(x)
}
