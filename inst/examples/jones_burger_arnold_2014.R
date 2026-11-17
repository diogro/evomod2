# Reproducing Jones, Burger & Arnold (2014), "Epistasis and natural selection
# shape the mutational architecture of complex traits", Nature Communications
# 5: 3709, doi:10.1038/ncomms4709.
#
# Model (their Methods and Table 1): 2 traits, n = 20 unlinked diploid loci,
# x = sum_i y_i + sum_{i<j} eps_ij(y_i, y_j) with a different random epistasis
# tensor for every pair of loci (mean 0, variance sigma_e^2), fixed in a run;
# mutation adds a bivariate normal vector (variances 0.05, correlation 0) to an
# allele's reference effect at rate mu = 0.0005 per locus; environmental
# variance 1; Gaussian stabilizing selection around a stationary optimum, with
# omega_11 = omega_22 = 49 and correlation r_omega; each locus starts with five
# equally frequent alleles with effects of standard deviation
# sqrt(0.05) / n; 5000 burn-in generations, then 5000 experimental generations
# in which G, the additive covariance and M are measured (20 times, averaged);
# 20 runs per setting.
#
# Life cycle: Wright-Fisher by default; with the environment variable
# EVOMOD2_LIFE_CYCLE=viability, their life cycle (viability_selection():
# females mating with two males, four offspring each, viability selection and
# regulation to K), saved as vs_table2_* and vs_table3_* files.
# Differences: the additive
# covariance is the least-squares one (additive_covariance()) rather than a
# half-sib analysis of variance.
#
# Experiments:
#   table2 <sigma_e^2> <r_omega> <replicate> <outdir>   N = 512 (their Table 2)
#   table3 <N> <r_omega> <replicate> <outdir>           sigma_e^2 = 1 (Table 3, Fig. 1)
#   fig2 <replicate> <outdir>     N = 4096, omega (r = 0.9 shape) rotated by replicate
#   fig3 <N> <replicate> <outdir> 10 loci, epistasis only within trait 1 (Fig. 3)
#   summarise <outdir> <summary.rds>

library(evomod2)

# ---- Settings (Jones, Burger & Arnold 2014, Table 1) -----------------------------

TRAITS <- 2
LOCI <- 20
N_CORE <- 512
ALPHA2 <- 0.05             # mutational variance per trait
R_MU <- 0                  # mutational correlation
MU <- 0.0005               # per locus (allele) per generation
OMEGA_VAR <- 49
ENV_SD <- 1
BURN_IN <- 5000
EXPERIMENT <- 5000
MEASUREMENTS <- 20         # measurements during the experimental phase
M_SAMPLES <- 5000          # single mutations per measurement
VIABILITY <- identical(Sys.getenv("EVOMOD2_LIFE_CYCLE"), "viability")
LIFE_CYCLE <- if (VIABILITY) viability_selection(mates = 2, offspring_per_mate = 2) else wright_fisher()
PREFIX <- if (VIABILITY) "vs_" else ""

omega_matrix <- function(r) OMEGA_VAR * matrix(c(1, r, r, 1), 2)
M_COV <- ALPHA2 * matrix(c(1, R_MU, R_MU, 1), 2)

# five equally frequent alleles per locus, effects with sd sqrt(alpha^2) / n
initial_genomes <- function(n, loci) {
  g <- array(0, c(TRAITS, 2 * loci, n))
  for (l in seq_len(loci)) {
    types <- matrix(stats::rnorm(5 * TRAITS, 0, sqrt(ALPHA2) / loci), TRAITS)
    g[, c(2 * l - 1, 2 * l), ] <- types[, sample(5, 2 * n, replace = TRUE)]
  }
  g
}

# ---- One run ------------------------------------------------------------------------

run <- function(n, sigma2, omega, loci = LOCI, eps = NULL, seed) {
  set.seed(seed)
  if (is.null(eps)) eps <- pairwise_epsilon(loci, TRAITS, sd = sqrt(sigma2))
  model <- evomod_model(multilinear(loci, TRAITS, eps), mutation(MU, M_COV), env_sd = ENV_SD,
                        life_cycle = LIFE_CYCLE)
  pop <- evolve(model, stabilizing(BURN_IN, omega = omega), n = n, record_every = BURN_IN,
                initial = initial_genomes(n, loci))
  stats <- lapply(seq_len(MEASUREMENTS), function(i) {
    pop <<- evolve(pop, stabilizing(EXPERIMENT / MEASUREMENTS, omega = omega))
    a <- additive_covariance(pop)
    mc <- mutational_covariance(pop, n = M_SAMPLES)
    list(G = a$G, A = a$additive, M = mc$per_mutation)
  })
  avg <- function(k) Reduce(`+`, lapply(stats, `[[`, k)) / MEASUREMENTS
  list(G = avg("G"), A = avg("A"), M = avg("M"), pop = pop, eps = eps)
}

cor12 <- function(m) m[1, 2] / sqrt(m[1, 1] * m[2, 2])
leading <- function(m) eigen(m, symmetric = TRUE)$vectors[, 1]
angle <- function(u, v) acos(min(1, abs(sum(u * v)))) * 180 / pi

matrix_row <- function(r) data.frame(
  VG11 = r$G[1, 1], VG22 = r$G[2, 2], rG = cor12(r$G),
  VA11 = r$A[1, 1], VA22 = r$A[2, 2], rA = cor12(r$A),
  M11 = r$M[1, 1], M22 = r$M[2, 2], rM = cor12(r$M))

seed_of <- function(...) sum(c(...) * c(1e6, 1e4, 1)[seq_along(c(...))])

run_table2 <- function(sigma2, r_omega, rep, outdir) {
  r <- run(N_CORE, sigma2, omega_matrix(r_omega), seed = seed_of(1, round(100 * sigma2) + round(10 * r_omega) * 1000, rep))
  saveRDS(data.frame(experiment = "table2", life_cycle = LIFE_CYCLE$type, N = N_CORE,
                     sigma2 = sigma2, r_omega = r_omega, replicate = rep, matrix_row(r)),
          file.path(outdir, sprintf("%stable2_%g_%g_%02d.rds", PREFIX, sigma2, r_omega, rep)))
}

run_table3 <- function(n, r_omega, rep, outdir) {
  r <- run(n, 1, omega_matrix(r_omega), seed = seed_of(2, n + round(100 * r_omega), rep))
  saveRDS(data.frame(experiment = "table3", life_cycle = LIFE_CYCLE$type, N = n, sigma2 = 1,
                     r_omega = r_omega, replicate = rep, matrix_row(r)),
          file.path(outdir, sprintf("%stable3_%d_%g_%02d.rds", PREFIX, n, r_omega, rep)))
}

run_fig2 <- function(rep, outdir) {
  phi <- (rep - 1) * pi / 20                              # orientation of omega
  R <- matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), 2)
  base <- diag(OMEGA_VAR * c(1.9, 0.1))                   # shape of r_omega = 0.9
  omega <- R %*% base %*% t(R)
  r <- run(4096, 1, omega, seed = seed_of(3, 0, rep))
  saveRDS(data.frame(experiment = "fig2", N = 4096, sigma2 = 1, r_omega = 0.9, replicate = rep,
                     omega_angle = phi * 180 / pi,
                     angle_G_omega = angle(leading(r$G), leading(omega)),
                     angle_M_omega = angle(leading(r$M), leading(omega)),
                     angle_M_G = angle(leading(r$M), leading(r$G)),
                     G_ratio = do.call(`/`, as.list(eigen(r$G)$values)),
                     M_ratio = do.call(`/`, as.list(eigen(r$M)$values)),
                     matrix_row(r)),
          file.path(outdir, sprintf("fig2_%02d.rds", rep)))
}

run_fig3 <- function(n, rep, outdir) {
  loci <- 10
  set.seed(seed_of(4, n, rep))
  eps <- pairwise_epsilon(loci, TRAITS, sd = 1, within_trait = TRUE, traits_with_epistasis = 1)
  r <- run(n, 1, omega_matrix(0), loci = loci, eps = eps, seed = seed_of(4, n, rep))
  g <- genomes(r$pop)
  # per locus: mean reference effect on trait 1 (sum of the two alleles) and the
  # mean of its epistatic coefficients with the other loci (trait 1)
  e <- r$eps[1, 1, 1, , ]
  e <- e + t(e)
  saveRDS(data.frame(experiment = "fig3", N = n, replicate = rep, locus = seq_len(loci),
                     effect = sapply(seq_len(loci), function(l) mean(g[1, 2 * l - 1, ] + g[1, 2 * l, ])),
                     abs_effect = sapply(seq_len(loci), function(l) mean(abs(g[1, 2 * l - 1, ] + g[1, 2 * l, ]))),
                     epistasis = rowSums(e) / (loci - 1)),
          file.path(outdir, sprintf("fig3_%d_%02d.rds", n, rep)))
}

summarise_runs <- function(outdir, file) {
  files <- list.files(outdir, "\\.rds$", full.names = TRUE)
  read <- function(prefix) do.call(rbind, lapply(grep(paste0("/", prefix, "_"), files, value = TRUE), readRDS))
  saveRDS(list(table2 = read("table2"), table3 = read("table3"), fig2 = read("fig2"),
               fig3 = read("fig3"), viability_table2 = read("vs_table2"),
               viability_table3 = read("vs_table3"),
               settings = list(traits = TRAITS, loci = LOCI, n = N_CORE, alpha2 = ALPHA2, mu = MU,
                               omega = OMEGA_VAR, env_sd = ENV_SD, burn_in = BURN_IN,
                               experiment = EXPERIMENT,
                               evomod2 = as.character(utils::packageVersion("evomod2")))),
          file, compress = "xz")
}

# ---- Command line -------------------------------------------------------------------

if (!interactive() && sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  cmd <- if (length(args)) args[1] else ""
  if (cmd != "summarise") dir.create(args[length(args)], showWarnings = FALSE, recursive = TRUE)
  switch(cmd,
    table2 = run_table2(as.numeric(args[2]), as.numeric(args[3]), as.integer(args[4]), args[5]),
    table3 = run_table3(as.integer(args[2]), as.numeric(args[3]), as.integer(args[4]), args[5]),
    fig2 = run_fig2(as.integer(args[2]), args[3]),
    fig3 = run_fig3(as.integer(args[2]), as.integer(args[3]), args[4]),
    summarise = summarise_runs(args[2], args[3]),
    stop("usage: jones_burger_arnold_2014.R table2|table3|fig2|fig3|summarise ..."))
}
