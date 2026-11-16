#' Plot recorded statistics
#'
#' Base-graphics plots of the recorded statistics through time: mean
#' genotypic values and optima (`"means"`), additive genetic variances
#' (`"variances"`), or the genetic correlation between two traits
#' (`"correlation"`). Regime boundaries are marked with vertical lines.
#'
#' @param x an `evomod_population`.
#' @param what which statistic to plot.
#' @param traits for `"correlation"`, the two traits to correlate.
#' @param ... passed to [graphics::matplot()].
#' @examples
#' model <- evomod_model(additive(loci = 20, traits = 2),
#'                       mutation(rate = 0.01, cov = 0.1), env_sd = 1, init_sd = 0.1)
#' pop <- evolve(model, list(stabilizing(50, omega = 10),
#'                           directional(50, delta_theta = c(0.02, 0), omega = 10)),
#'               n = 100, record_every = 5, seed = 1)
#' plot(pop)
#' plot(pop, "variances")
#' @export
plot.evomod_population <- function(x, what = c("means", "variances", "correlation"),
                                   traits = c(1, 2), ...) {
  what <- match.arg(what)
  r <- x$records
  gen <- r$generation
  p <- x$model$traits
  cols <- grDevices::hcl.colors(p, "Dark 3")
  if (what == "means") {
    graphics::matplot(gen, t(r$mean_x), type = "l", lty = 1, col = cols,
                      xlab = "Generation", ylab = "Mean genotypic value", ...)
    graphics::matlines(gen, t(r$theta), lty = 2, col = cols)
    graphics::legend("topleft", c(paste0("trait", seq_len(p)), "optimum"),
                     col = c(cols, "grey40"), lty = c(rep(1, p), 2), bty = "n")
  } else if (what == "variances") {
    g <- matrix(apply(r$G, 3, diag), nrow = p)
    graphics::matplot(gen, t(g), type = "l", lty = 1, col = cols,
                      xlab = "Generation", ylab = "Additive genetic variance", ...)
    graphics::legend("topleft", paste0("trait", seq_len(p)), col = cols, lty = 1, bty = "n")
  } else {
    i <- traits[1]; j <- traits[2]
    rho <- r$G[i, j, ] / sqrt(r$G[i, i, ] * r$G[j, j, ])
    graphics::plot(gen, rho, type = "l", ylim = c(-1, 1), xlab = "Generation",
                   ylab = paste("Genetic correlation, traits", i, "and", j), ...)
  }
  starts <- vapply(x$regimes, `[[`, 1, "start")
  graphics::abline(v = starts[-1], col = "grey80")
  invisible(x)
}
