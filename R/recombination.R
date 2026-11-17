#' Recombination
#'
#' How the alleles a parent passes on are put together:
#'
#' * `free_recombination()` (the default): every locus is inherited
#'   independently.
#' * `genetic_map()`: loci are placed in order, evenly spaced, on `chromosomes`
#'   chromosomes of equal length (`length`, in Morgans; loci `1..k` on the first
#'   chromosome, and so on). Recombination fractions between neighbouring loci
#'   follow Haldane's map function, \eqn{r = (1 - e^{-2d}) / 2}{(1 - exp(-2 d)) / 2} for a
#'   distance \eqn{d}{d};
#'   loci on different chromosomes recombine with probability 1/2.
#' * `linkage()`: explicit recombination fractions between neighbouring loci,
#'   a vector of length `loci - 1` (0.5 for unlinked loci).
#'
#' Columns `2l - 1` and `2l` of a genome are the paternal and maternal
#' haplotypes at locus `l`; a gamete starts on a random haplotype and switches
#' at each crossover.
#'
#' @param chromosomes number of chromosomes.
#' @param length length of each chromosome in Morgans (1 Morgan = 100 cM).
#' @param rates recombination fractions between loci `l` and `l + 1`.
#' @return An `evomod_recombination` object, resolved for the number of loci
#'   by [evomod_model()].
#' @examples
#' genetic_map(chromosomes = 2, length = 0.5)
#' linkage(c(0.01, 0.5, 0.01))
#' @name recombination
NULL

#' @rdname recombination
#' @export
free_recombination <- function() structure(list(type = "free"), class = "evomod_recombination")

#' @rdname recombination
#' @export
genetic_map <- function(chromosomes = 1, length = 1) {
  check_count(chromosomes, "chromosomes")
  if (!is.numeric(length) || length(length) != 1 || length <= 0)
    stop("`length` must be a positive number of Morgans", call. = FALSE)
  structure(list(type = "genetic_map", chromosomes = as.integer(chromosomes), length = length),
            class = "evomod_recombination")
}

#' @rdname recombination
#' @export
linkage <- function(rates) {
  if (!is.numeric(rates) || any(rates < 0 | rates > 0.5))
    stop("`rates` must be recombination fractions between 0 and 0.5", call. = FALSE)
  structure(list(type = "linkage", rates = as.numeric(rates)), class = "evomod_recombination")
}

haldane <- function(d) (1 - exp(-2 * d)) / 2

# Recombination as the C++ code expects it: "free", or "map" with the
# recombination fraction between each pair of neighbouring loci.
resolve_recombination <- function(recombination, loci) {
  if (!inherits(recombination, "evomod_recombination"))
    stop("`recombination` must be created with free_recombination(), genetic_map() or linkage()",
         call. = FALSE)
  if (recombination$type == "free") return(recombination)
  rates <- switch(recombination$type,
    linkage = {
      if (length(recombination$rates) != loci - 1)
        stop("`rates` must have one value per pair of neighbouring loci (", loci - 1, ")",
             call. = FALSE)
      recombination$rates
    },
    genetic_map = {
      k <- recombination$chromosomes
      if (k > loci) stop("more chromosomes than loci", call. = FALSE)
      chrom <- sort(rep_len(seq_len(k), loci))
      per_chrom <- tabulate(chrom, k)
      rates <- numeric(loci - 1)
      for (l in seq_len(loci - 1)) {
        rates[l] <- if (chrom[l] != chrom[l + 1]) 0.5 else
          haldane(recombination$length / (per_chrom[chrom[l]] - 1))
      }
      recombination$chromosome <- chrom
      rates
    })
  recombination$rates <- rates
  recombination$map_type <- recombination$type
  recombination$type <- "map"
  recombination
}

#' @export
print.evomod_recombination <- function(x, ...) {
  cat("<evomod recombination>",
      switch(x$type, free = "free recombination",
             genetic_map = paste(x$chromosomes, "chromosome(s) of", x$length, "Morgans"),
             linkage = paste("recombination fractions for", length(x$rates) + 1, "loci"),
             map = paste("map of", length(x$rates) + 1, "loci")), "\n")
  invisible(x)
}
