# Reproducing Melo & Marroig (2015), "Directional selection can drive the
# evolution of modularity in complex traits", PNAS 112(2): 470-475,
# doi:10.1073/pnas.1322632112.
#
# Model (Methods of the paper; C code at github.com/lem-usp/evomod):
#   10 traits in two potential modules (traits 1-5 and 6-10), 500 loci with a
#   binary ontogenetic matrix B, N = 5000, mu = 5e-4 per allele,
#   mu_B = 1e-4 per entry of B, mutational effects with standard deviation
#   0.02 (called a variance in the paper; the C code used it as a standard
#   deviation), environmental standard deviation 0.8 (likewise), Gaussian
#   fitness with V_omega = 10. The published simulations set a mutated B entry
#   to a random value (b_mutation = "random").
#
# Protocol:
#   burnin       10,000 generations of drift, then 10,000 of uncorrelated
#                stabilizing selection, from a random population; one
#                population is the starting point of every experiment
#   divergent    10,000 generations of divergent directional selection
#                (traits 1-5 up, 6-10 down) on a correlated fitness surface
#                (r_omega = 0.8 within modules), Fig. 2
#   corridor     10,000 generations of directional selection on module 1
#                (traits 1-5), module 2 under stabilizing selection, Fig. 3
#   maintenance  10,000 generations of drift, uncorrelated or correlated
#                stabilizing selection, starting from a population after
#                strong divergent selection (rate 0.004), Fig. 4
#   summarise    collect the statistics of every run in one small file
#
# Each run takes about 40 minutes and 2 GB of memory; run replicates in
# parallel, one process each:
#   Rscript melo_marroig_2015.R burnin <outdir>
#   Rscript melo_marroig_2015.R divergent <rate> <replicate> <outdir>
#   Rscript melo_marroig_2015.R corridor <rate> <replicate> <outdir>
#   Rscript melo_marroig_2015.R maintenance <drift|uncorrelated|correlated> <replicate> <outdir>
#   Rscript melo_marroig_2015.R summarise <outdir> <summary.rds>
#
# The vignette "Reproducing Melo & Marroig 2015" plots the summary shipped
# with the package (system.file("extdata", "melo_marroig_2015.rds",
# package = "evomod2")).

library(evomod2)

# ---- Settings (Melo & Marroig 2015) --------------------------------------------

TRAITS <- 10
LOCI <- 500
N <- 5000
MU <- 5e-4
MU_B <- 1e-4
SIGMA <- 0.02
V_E <- 0.8
V_OMEGA <- 10
R_OMEGA <- 0.8
MODULES <- list(1:5, 6:10)
GENERATIONS <- 10000
RECORD_EVERY <- 100
RATES <- seq(0.0004, 0.004, by = 0.0004)   # peak movement per generation
MAINTENANCE_SEED_RATE <- 0.004

model <- evomod_model(
  binary_b(LOCI, TRAITS, init_b = 0.25, init_y = "mutations"),
  binary_b_mutation(rate = MU, sd = SIGMA, rate_b = MU_B, b_mutation = "random"),
  env_sd = V_E, init_sd = SIGMA)

omega_uncorrelated <- diag(V_OMEGA, TRAITS)
omega_correlated <- V_OMEGA * (kronecker(diag(2), matrix(R_OMEGA, 5, 5)) + diag(1 - R_OMEGA, TRAITS))

# ---- Statistics of a run ----------------------------------------------------------

module_1 <- rep(c(1, 0), each = 5)

matrix_stats <- function(m) {
  a <- avg_corr(m, modules = MODULES)
  c(avg_ratio = a$avg_ratio, within_1 = a$within[1], within_2 = a$within[2],
    between = a$between,
    setNames(eigen_variance(m)[1:4], paste0("eigen_", 1:4)),
    # NA if the matrix is singular (a trait without variation)
    autonomy_module = tryCatch(autonomy(m, module_1), error = function(e) NA),
    autonomy_mean = tryCatch(mean_autonomy(m, 1000), error = function(e) NA))
}

# Statistics of P (as in the paper) and G at every recorded generation after
# `from`, and the final P and G matrices.
run_summary <- function(pop, from, info) {
  P <- p_matrices(pop)
  G <- g_matrices(pop)
  gens <- as.integer(names(P))
  keep <- gens >= from
  trace <- do.call(rbind, lapply(which(keep), function(i)
    data.frame(info, generation = gens[i] - from, matrix = c("P", "G"),
               rbind(matrix_stats(P[[i]]), matrix_stats(G[[i]])))))
  list(info = info, trace = trace, final_P = P[[length(P)]], final_G = G[[length(G)]],
       h2 = variances(pop)[variances(pop)$generation == max(gens), ])
}

path <- function(outdir, ...) file.path(outdir, paste0(...))

# ---- Experiments ----------------------------------------------------------------

run_burnin <- function(outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  pop <- evolve(model, n = N, record_every = RECORD_EVERY, seed = 2015,
                regimes = list(drift(GENERATIONS),
                               stabilizing(GENERATIONS, omega = omega_uncorrelated)))
  saveRDS(pop, path(outdir, "burnin_population.rds"))
  saveRDS(run_summary(pop, 0, data.frame(experiment = "burnin", rate = 0, replicate = 1)),
          path(outdir, "burnin.rds"))
}

start_population <- function(outdir) readRDS(path(outdir, "burnin_population.rds"))

run_divergent <- function(rate, replicate, outdir) {
  pop <- start_population(outdir)
  start <- pop$generation
  pop <- evolve(pop, directional(GENERATIONS, delta_theta = rate * rep(c(1, -1), each = 5),
                                 omega = omega_correlated),
                seed = 1e5 + round(rate * 1e4) * 100 + replicate)
  if (isTRUE(all.equal(rate, MAINTENANCE_SEED_RATE)) && replicate == 1)
    saveRDS(pop, path(outdir, "maintenance_seed_population.rds"))
  saveRDS(run_summary(pop, start, data.frame(experiment = "divergent", rate = rate, replicate = replicate)),
          path(outdir, sprintf("divergent_%.4f_%02d.rds", rate, replicate)))
}

run_corridor <- function(rate, replicate, outdir) {
  pop <- start_population(outdir)
  start <- pop$generation
  pop <- evolve(pop, directional(GENERATIONS, delta_theta = rate * rep(c(1, 0), each = 5),
                                 omega = omega_correlated),
                seed = 2e5 + round(rate * 1e4) * 100 + replicate)
  saveRDS(run_summary(pop, start, data.frame(experiment = "corridor", rate = rate, replicate = replicate)),
          path(outdir, sprintf("corridor_%.4f_%02d.rds", rate, replicate)))
}

run_maintenance <- function(regime, replicate, outdir) {
  pop <- readRDS(path(outdir, "maintenance_seed_population.rds"))
  start <- pop$generation
  reg <- switch(regime,
    drift = drift(GENERATIONS),
    uncorrelated = stabilizing(GENERATIONS, omega = omega_uncorrelated),
    correlated = stabilizing(GENERATIONS, omega = omega_correlated),
    stop("regime must be drift, uncorrelated or correlated"))
  pop <- evolve(pop, reg, seed = 3e5 + match(regime, c("drift", "uncorrelated", "correlated")) * 100 + replicate)
  saveRDS(run_summary(pop, start, data.frame(experiment = "maintenance", rate = 0, replicate = replicate,
                                             regime = regime)),
          path(outdir, sprintf("maintenance_%s_%02d.rds", regime, replicate)))
}

summarise_runs <- function(outdir, file) {
  files <- setdiff(list.files(outdir, "\\.rds$", full.names = TRUE),
                   path(outdir, c("burnin_population.rds", "maintenance_seed_population.rds")))
  runs <- lapply(files, readRDS)
  fill <- function(d) { if (is.null(d$regime)) d$regime <- NA_character_; d }
  trace <- do.call(rbind, lapply(runs, function(r) fill(r$trace)))
  # keep full traces for the burn-in and maintenance runs; the final
  # generation for the rest
  last <- ave(trace$generation, trace$experiment, trace$rate, trace$replicate, trace$regime,
              FUN = max)
  trace <- trace[trace$experiment %in% c("burnin", "maintenance") | trace$generation == last, ]
  burnin <- runs[[which(vapply(runs, function(r) r$info$experiment == "burnin", TRUE))]]
  summary <- list(
    statistics = trace,
    burnin_h2 = burnin$h2,
    final_P = lapply(runs, function(r) list(info = fill(r$info), P = r$final_P, G = r$final_G)),
    settings = list(traits = TRAITS, loci = LOCI, n = N, mu = MU, mu_b = MU_B, sigma = SIGMA,
                    v_e = V_E, v_omega = V_OMEGA, r_omega = R_OMEGA, generations = GENERATIONS,
                    record_every = RECORD_EVERY, rates = RATES,
                    evomod2 = as.character(utils::packageVersion("evomod2"))))
  saveRDS(summary, file, compress = "xz")
}

# ---- Command line -----------------------------------------------------------------

if (!interactive() && sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  cmd <- if (length(args)) args[1] else ""
  switch(cmd,
    burnin = run_burnin(args[2]),
    divergent = run_divergent(as.numeric(args[2]), as.integer(args[3]), args[4]),
    corridor = run_corridor(as.numeric(args[2]), as.integer(args[3]), args[4]),
    maintenance = run_maintenance(args[2], as.integer(args[3]), args[4]),
    summarise = summarise_runs(args[2], args[3]),
    stop("usage: melo_marroig_2015.R burnin|divergent|corridor|maintenance|summarise ..."))
}
