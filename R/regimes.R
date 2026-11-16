#' Selection regimes
#'
#' A simulation runs through a list of regimes:
#'
#' * `drift()`: no selection.
#' * `stabilizing()`: Gaussian fitness surface with optimum `theta` and
#'   covariance `omega`, `w(z) = exp(-(z - theta)' omega^-1 (z - theta) / 2)`.
#' * `directional()`: as stabilizing, but the optimum moves by `delta_theta`
#'   every generation (before selection).
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
#' @return An `evomod_regime` object.
#' @examples
#' regimes <- list(drift(100), stabilizing(1000, omega = 20),
#'                 directional(1000, delta_theta = c(0.005, 0), omega = 20))
#' @name regimes
NULL

new_regime <- function(type, generations, omega = NULL, theta = NULL, delta_theta = NULL) {
  check_count(generations, "generations")
  structure(list(type = type, generations = as.integer(generations), omega = omega,
                 theta = theta, delta_theta = delta_theta),
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

#' @export
print.evomod_regime <- function(x, ...) {
  cat("<evomod regime>", x$type, "for", x$generations, "generations\n")
  invisible(x)
}

# Regimes as the C++ code expects them, with the optimum at the start of each
# regime chained from the previous one.
resolve_regimes <- function(regimes, p, theta) {
  if (inherits(regimes, "evomod_regime")) regimes <- list(regimes)
  if (!length(regimes) || !all(vapply(regimes, inherits, TRUE, "evomod_regime")))
    stop("`regimes` must be a list of drift(), stabilizing() or directional() regimes",
         call. = FALSE)
  lapply(regimes, function(r) {
    omega <- if (is.null(r$omega)) diag(p) else r$omega
    if (length(omega) == 1 && !is.matrix(omega)) omega <- diag(omega, p)
    omega <- as.matrix(omega)
    if (!identical(dim(omega), c(p, p)) || !isSymmetric(unname(omega)))
      stop("`omega` must be a symmetric ", p, " x ", p, " matrix", call. = FALSE)
    if (!is.null(r$theta)) theta <<- check_vector(r$theta, p, "theta")
    delta <- if (is.null(r$delta_theta)) rep(0, p) else check_vector(r$delta_theta, p, "delta_theta")
    out <- list(generations = r$generations, selective = r$type != "drift",
                theta = theta, delta_theta = delta, omega = omega, type = r$type)
    theta <<- theta + r$generations * delta
    out
  })
}

check_vector <- function(x, p, name) {
  if (!is.numeric(x) || length(x) != p)
    stop("`", name, "` must be a numeric vector of length ", p, call. = FALSE)
  as.numeric(x)
}
