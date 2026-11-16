# Corridor selection with and without linkage: pleiotropy vs. linkage
# disequilibrium as sources of genetic correlation.
#
# Same setting as corridor_selection.R (four traits, modules 1-2 and 3-4;
# drift, stabilizing selection, then corridor selection on module 1), for the
# binary B and additive architectures, each with free recombination and with
# tight linkage (all loci on one short chromosome). G is split into its genic
# part (pleiotropy: allele effects on several traits) and its linkage
# disequilibrium part (associations between loci built by selection), see
# ?g_decomposition. Under free recombination LD decays quickly; under tight
# linkage it persists, so an additive genome can build correlations without
# pleiotropy.
#
# Each part of G is expressed on the correlation scale of the total G
# (C_ij / sqrt(G_ii G_jj)), so genic + LD = total correlation.
#
# Usage (from a shell):
#   Rscript corridor_linkage.R                          # everything, sequentially
#   Rscript corridor_linkage.R simulate <condition> <replicate> <outdir>
#   Rscript corridor_linkage.R plot <outdir>
# Conditions: binary_b_free, binary_b_linked, additive_free, additive_linked.

library(evomod2)

# ---- Settings ------------------------------------------------------------------

TRAITS <- 4
LOCI <- 100
N <- 2000
GENERATIONS <- c(drift = 1000, stabilizing = 5000, corridor = 5000)
RECORD_EVERY <- 50
REPLICATES <- 20

ENV_SD <- 0.8
MU <- 0.0005
SIGMA <- 0.2
MU_B <- 0.0001

OMEGA <- diag(10, TRAITS)
SPEED <- 0.002
DELTA_THETA <- c(SPEED, SPEED, 0, 0)
MODULES <- list(1:2, 3:4)

TIGHT <- genetic_map(chromosomes = 1, length = 0.5)   # ~0.005 between neighbouring loci

model <- function(arch, recombination) {
  switch(arch,
    binary_b = evomod_model(binary_b(LOCI, TRAITS),
                            binary_b_mutation(rate = MU, sd = SIGMA, rate_b = MU_B),
                            env_sd = ENV_SD, init_sd = SIGMA, recombination = recombination),
    additive = evomod_model(additive(LOCI, TRAITS), mutation(rate = MU, cov = SIGMA),
                            env_sd = ENV_SD, init_sd = SIGMA, recombination = recombination))
}

conditions <- list(
  binary_b_free   = list(arch = "binary_b", recombination = "free"),
  binary_b_linked = list(arch = "binary_b", recombination = "linked"),
  additive_free   = list(arch = "additive", recombination = "free"),
  additive_linked = list(arch = "additive", recombination = "linked"))

regimes <- list(drift(GENERATIONS[["drift"]]),
                stabilizing(GENERATIONS[["stabilizing"]], omega = OMEGA),
                directional(GENERATIONS[["corridor"]], delta_theta = DELTA_THETA, omega = OMEGA))

# ---- Simulation ----------------------------------------------------------------

# average within-module and between-module values of C on the correlation
# scale of G
scaled_avg <- function(C, G) {
  R <- C / sqrt(outer(diag(G), diag(G)))
  within <- vapply(MODULES, function(m) mean(R[m, m][lower.tri(R[m, m])]), 1)
  same <- matrix(FALSE, TRAITS, TRAITS)
  for (m in MODULES) same[m, m] <- TRUE
  c(module_1 = within[1], module_2 = within[2], between = mean(R[!same]))
}

simulate_replicate <- function(condition, replicate) {
  cond <- conditions[[condition]]
  rec <- if (cond$recombination == "free") free_recombination() else TIGHT
  pop <- evolve(model(cond$arch, rec), regimes, n = N, record_every = RECORD_EVERY,
                seed = replicate)
  G <- g_matrices(pop)
  parts <- list(total = G, genic = g_matrices(pop, "genic"), ld = g_matrices(pop, "ld"))
  corr <- do.call(rbind, lapply(names(G), function(gen) {
    do.call(rbind, lapply(names(parts), function(part) {
      v <- scaled_avg(parts[[part]][[gen]], G[[gen]])
      data.frame(condition = condition, architecture = cond$arch,
                 recombination = cond$recombination, replicate = replicate,
                 generation = as.integer(gen), part = part,
                 correlation = names(v), value = unname(v))
    }))
  }))
  list(correlations = corr, final = g_decomposition(pop))
}

save_replicate <- function(condition, replicate, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(simulate_replicate(condition, replicate),
          file.path(outdir, sprintf("%s_%03d.rds", condition, replicate)))
}

# ---- Plots ---------------------------------------------------------------------

plot_results <- function(outdir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("install ggplot2 for the plots")
  library(ggplot2)
  files <- list.files(outdir, pattern = "\\.rds$", full.names = TRUE)
  corr <- do.call(rbind, lapply(files, function(f) readRDS(f)$correlations))
  corr$correlation <- factor(corr$correlation, c("module_1", "module_2", "between"),
                             c("Within module 1\n(selected)", "Within module 2", "Between modules"))
  corr$architecture <- factor(corr$architecture, c("binary_b", "additive"),
                              c("Binary B matrix", "Additive"))
  corr$recombination <- factor(corr$recombination, c("free", "linked"),
                               c("Free recombination", "Tight linkage"))
  boundaries <- cumsum(GENERATIONS)[-length(GENERATIONS)]
  colours <- c("#D55E00", "#0072B2", "#7F7F7F")
  n_rep <- length(unique(corr$replicate))

  # total correlations through time
  total <- corr[corr$part == "total", ]
  mean_total <- aggregate(value ~ architecture + recombination + correlation + generation, total, mean)
  trace <- ggplot(total, aes(generation, value, colour = correlation)) +
    geom_vline(xintercept = boundaries, colour = "grey85") +
    geom_line(aes(group = interaction(replicate, correlation)), alpha = 0.12, linewidth = 0.3) +
    geom_line(data = mean_total, linewidth = 0.9) +
    facet_grid(architecture ~ recombination) +
    scale_colour_manual(values = colours, name = NULL) +
    labs(x = "Generation", y = "Average genetic correlation",
         title = "Corridor selection on module 1, with and without linkage",
         subtitle = sprintf("%d replicates per panel; drift, stabilizing and corridor phases separated by grey lines", n_rep)) +
    theme_minimal() + theme(legend.position = "bottom")

  # genic and LD parts through time (means over replicates)
  parts <- corr[corr$part != "total", ]
  parts$part <- factor(parts$part, c("genic", "ld"), c("Genic (pleiotropy)", "Linkage disequilibrium"))
  mean_parts <- aggregate(value ~ architecture + recombination + correlation + part + generation,
                          parts, mean)
  parts_trace <- ggplot(mean_parts, aes(generation, value, colour = correlation, linetype = part)) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = boundaries, colour = "grey85") +
    geom_line(linewidth = 0.7) +
    facet_grid(architecture ~ recombination) +
    scale_colour_manual(values = colours, name = NULL) +
    scale_linetype_manual(values = c("solid", "22"), name = NULL) +
    labs(x = "Generation", y = "Contribution to average genetic correlation",
         title = "Genic and linkage disequilibrium parts of the genetic correlations",
         subtitle = sprintf("Means of %d replicates; the two parts add up to the total correlation", n_rep)) +
    theme_minimal() + theme(legend.position = "bottom", legend.box = "vertical")

  # final generation: stacked genic + LD contributions
  final <- parts[parts$generation == max(parts$generation), ]
  final_mean <- aggregate(value ~ architecture + recombination + correlation + part, final, mean)
  final_plot <- ggplot(final_mean, aes(correlation, value, fill = part)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_col(width = 0.6) +
    facet_grid(architecture ~ recombination) +
    scale_fill_manual(values = c("#009E73", "#CC79A7"), name = NULL) +
    labs(x = NULL, y = "Average genetic correlation",
         title = "Final G matrices: genic and linkage disequilibrium parts",
         subtitle = sprintf("Generation %d, means of %d replicates", max(final$generation), n_rep)) +
    theme_minimal() + theme(legend.position = "bottom")

  ggsave(file.path(outdir, "linkage_correlation_traces.png"), trace, width = 9, height = 7, dpi = 150)
  ggsave(file.path(outdir, "linkage_parts_traces.png"), parts_trace, width = 9, height = 7, dpi = 150)
  ggsave(file.path(outdir, "linkage_final_parts.png"), final_plot, width = 9, height = 7, dpi = 150)
  summary <- reshape(final_mean, idvar = c("architecture", "recombination", "correlation"),
                     timevar = "part", direction = "wide")
  names(summary) <- sub("value.", "", names(summary), fixed = TRUE)
  summary$total <- summary[["Genic (pleiotropy)"]] + summary[["Linkage disequilibrium"]]
  summary$correlation <- gsub("\n", " ", summary$correlation)
  print(summary, digits = 2, row.names = FALSE)
  invisible(list(trace = trace, parts = parts_trace, final = final_plot, summary = summary))
}

# ---- Command line ---------------------------------------------------------------

if (!interactive() && sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    outdir <- "corridor_linkage_results"
    for (cond in names(conditions)) for (r in seq_len(REPLICATES)) save_replicate(cond, r, outdir)
    plot_results(outdir)
  } else if (args[1] == "simulate") {
    save_replicate(args[2], as.integer(args[3]), args[4])
  } else if (args[1] == "plot") {
    plot_results(args[2])
  } else {
    stop("usage: corridor_linkage.R [simulate <condition> <replicate> <outdir> | plot <outdir>]")
  }
}
