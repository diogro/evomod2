# Corridor selection on a two-module, four-trait system:
# binary B matrix vs. additive genetic architecture.
#
# Traits 1-2 form module 1 and traits 3-4 module 2. After drift and
# stabilizing selection, the optimum of module 1 moves in a straight line
# (corridor selection) while module 2 stays under stabilizing selection
# (Melo & Marroig 2015). The two architectures differ in how pleiotropy can
# evolve:
#
#   binary_b(): each allele affects the traits marked in its column of B, and B
#               evolves, so the pattern of pleiotropy itself can change.
#   additive(): each allele has an effect vector with independent mutations
#               in each trait, so genetic correlations come only from
#               linkage disequilibrium built by selection.
#
# The script records G matrices through time and plots the average genetic
# correlation within each module and between modules (trace plots), for
# replicate populations of each architecture.
#
# Usage (from a shell):
#   Rscript corridor_selection.R                     # everything, sequentially
#   Rscript corridor_selection.R simulate <arch> <replicate> <outdir>
#                                                    # one replicate (arch: binary_b or additive)
#   Rscript corridor_selection.R plot <outdir>       # figures from saved replicates
#
# Or from R: source(system.file("examples", "corridor_selection.R", package = "evomod2"))
#
# The full settings take a few minutes per replicate. To try it quickly, set
# N, GENERATIONS and REPLICATES smaller.

library(evomod2)

# ---- Settings ------------------------------------------------------------------

TRAITS <- 4
LOCI <- 100
N <- 2000                       # population size
GENERATIONS <- c(drift = 1000, stabilizing = 5000, corridor = 5000)
RECORD_EVERY <- 50
REPLICATES <- 20

ENV_SD <- 0.8                   # environmental standard deviation
MU <- 0.0005                    # mutation rate of allelic effects, per allele
SIGMA <- 0.2                    # standard deviation of mutational effects
MU_B <- 0.0001                  # binary B: mutation rate per entry of B

OMEGA <- diag(10, TRAITS)       # fitness surface
SPEED <- 0.002                  # movement of the optimum per generation
DELTA_THETA <- c(SPEED, SPEED, 0, 0)    # corridor: only module 1 moves
MODULES <- list(1:2, 3:4)

architectures <- list(
  binary_b = evomod_model(binary_b(LOCI, TRAITS),
                          binary_b_mutation(rate = MU, sd = SIGMA, rate_b = MU_B),
                          env_sd = ENV_SD, init_sd = SIGMA),
  additive = evomod_model(additive(LOCI, TRAITS),
                          mutation(rate = MU, cov = SIGMA),
                          env_sd = ENV_SD, init_sd = SIGMA)
)

regimes <- list(drift(GENERATIONS[["drift"]]),
                stabilizing(GENERATIONS[["stabilizing"]], omega = OMEGA),
                directional(GENERATIONS[["corridor"]], delta_theta = DELTA_THETA, omega = OMEGA))

# ---- Simulation ----------------------------------------------------------------

# One replicate: correlations within and between modules for every recorded G
simulate_replicate <- function(arch, replicate) {
  pop <- evolve(architectures[[arch]], regimes, n = N, record_every = RECORD_EVERY,
                seed = replicate)
  G <- g_matrices(pop)
  corr <- do.call(rbind, lapply(names(G), function(gen) {
    a <- avg_corr(G[[gen]], modules = MODULES)
    data.frame(architecture = arch, replicate = replicate, generation = as.integer(gen),
               module_1 = a$within[1], module_2 = a$within[2], between = a$between)
  }))
  list(correlations = corr, final_G = G[[length(G)]], means = trait_means(pop))
}

save_replicate <- function(arch, replicate, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  res <- simulate_replicate(arch, replicate)
  saveRDS(res, file.path(outdir, sprintf("%s_%03d.rds", arch, replicate)))
}

# ---- Plots ---------------------------------------------------------------------

plot_results <- function(outdir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("install ggplot2 for the plots")
  library(ggplot2)
  files <- list.files(outdir, pattern = "\\.rds$", full.names = TRUE)
  res <- lapply(files, readRDS)
  corr <- do.call(rbind, lapply(res, `[[`, "correlations"))
  long <- do.call(rbind, lapply(c("module_1", "module_2", "between"), function(k)
    data.frame(corr[, c("architecture", "replicate", "generation")],
               correlation = k, value = corr[[k]])))
  long$correlation <- factor(long$correlation, c("module_1", "module_2", "between"),
                             c("Within module 1 (selected)", "Within module 2",
                               "Between modules"))
  long$architecture <- factor(long$architecture, c("binary_b", "additive"),
                              c("Binary B matrix", "Additive"))
  boundaries <- cumsum(GENERATIONS)[-length(GENERATIONS)]
  colours <- c("#D55E00", "#0072B2", "#7F7F7F")

  # mean over replicates through time
  mean_corr <- aggregate(value ~ architecture + correlation + generation, long, mean)
  trace <- ggplot(long, aes(generation, value, colour = correlation)) +
    geom_vline(xintercept = boundaries, colour = "grey85") +
    geom_line(aes(group = interaction(replicate, correlation)), alpha = 0.15, linewidth = 0.3) +
    geom_line(data = mean_corr, linewidth = 0.9) +
    facet_wrap(~ architecture) +
    scale_colour_manual(values = colours, name = NULL) +
    labs(x = "Generation", y = "Average genetic correlation",
         title = "Corridor selection on module 1",
         subtitle = sprintf("%d replicates; drift, stabilizing and corridor phases separated by grey lines",
                            length(unique(long$replicate)))) +
    theme_minimal() + theme(legend.position = "bottom")

  # final G: average correlations at the last generation
  final <- long[long$generation == max(long$generation), ]
  final_plot <- ggplot(final, aes(correlation, value, colour = correlation)) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.12, alpha = 0.5, size = 1) +
    facet_wrap(~ architecture) +
    scale_colour_manual(values = colours, guide = "none") +
    scale_x_discrete(labels = function(x) sub(" (", "\n(", sub("Within ", "Within\n", x), fixed = TRUE)) +
    labs(x = NULL, y = "Average genetic correlation",
         title = "Final G matrices", subtitle = sprintf("Generation %d", max(long$generation))) +
    theme_minimal()

  # mean final G correlation matrix of each architecture
  final_G <- lapply(split(res, vapply(res, function(r) r$correlations$architecture[1], "")),
                    function(rs) Reduce(`+`, lapply(rs, function(r) cov2cor(r$final_G))) / length(rs))

  ggsave(file.path(outdir, "correlation_traces.png"), trace, width = 9, height = 4.5, dpi = 150)
  ggsave(file.path(outdir, "final_correlations.png"), final_plot, width = 7, height = 4.5, dpi = 150)
  print(lapply(final_G, round, 2))
  invisible(list(trace = trace, final = final_plot, final_G = final_G))
}

# ---- Command line ---------------------------------------------------------------

if (!interactive() && sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    outdir <- "corridor_results"
    for (arch in names(architectures)) for (r in seq_len(REPLICATES))
      save_replicate(arch, r, outdir)
    plot_results(outdir)
  } else if (args[1] == "simulate") {
    save_replicate(args[2], as.integer(args[3]), args[4])
  } else if (args[1] == "plot") {
    plot_results(args[2])
  } else {
    stop("usage: corridor_selection.R [simulate <arch> <replicate> <outdir> | plot <outdir>]")
  }
}
