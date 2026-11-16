#' Matrix statistics: autonomy and eigenvalue distribution
#'
#' * `autonomy()`: autonomy of a covariance matrix in the direction(s)
#'   `beta` (Hansen & Houle 2008), `a(beta) = 1 / ((beta' G^-1 beta) (beta' G beta))`
#'   for unit vectors `beta`: the fraction of the variation in direction
#'   `beta` that is independent of other directions (1 when `beta` is an
#'   eigenvector).
#' * `mean_autonomy()`: mean autonomy over `n` random directions.
#' * `eigen_variance()`: proportion of the total variance on each principal
#'   component.
#'
#' @param m a covariance matrix (e.g. a G or P matrix).
#' @param beta a direction (vector) or several directions (columns of a
#'   matrix); normalised to unit length.
#' @param n number of random directions.
#' @return `autonomy()`: one value per direction. `mean_autonomy()`: a number.
#'   `eigen_variance()`: a vector of proportions in decreasing order.
#' @examples
#' G <- matrix(c(1, 0.8, 0.1, 0.8, 1, 0.1, 0.1, 0.1, 1), 3)
#' autonomy(G, c(1, 1, 0))
#' mean_autonomy(G)
#' eigen_variance(G)
#' @name matrix_statistics
NULL

#' @rdname matrix_statistics
#' @export
autonomy <- function(m, beta) {
  m <- as.matrix(m)
  beta <- as.matrix(beta)
  beta <- sweep(beta, 2, sqrt(colSums(beta^2)), "/")
  inv <- solve(m)
  1 / (colSums(beta * (inv %*% beta)) * colSums(beta * (m %*% beta)))
}

#' @rdname matrix_statistics
#' @export
mean_autonomy <- function(m, n = 1000) {
  beta <- matrix(stats::rnorm(nrow(m) * n), nrow(m))
  mean(autonomy(m, beta))
}

#' @rdname matrix_statistics
#' @export
eigen_variance <- function(m) {
  ev <- eigen(as.matrix(m), symmetric = TRUE, only.values = TRUE)$values
  ev / sum(ev)
}
