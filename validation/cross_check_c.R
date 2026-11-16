# Cross-check the binary B architecture against the original C simulation
# (c-gsl in https://github.com/lem-usp/evomod), the model of Melo & Marroig (2015).
#
# Runs the same parameters in both and compares replicate means of the final G
# matrix and mean genotypic values. The implementations use different random
# number generators, so results agree in distribution, not exactly.
#
# Usage: Rscript validation/cross_check_c.R <path to evomod/c-gsl> [reps] [cores]
# Needs GSL (gsl-config) and a C compiler to build the C code. The C code
# flips mutated B entries since November 2026 (evomod 77ace1c); for older
# versions use b_mutation = "random" below.

args <- commandArgs(trailingOnly = TRUE)
c_dir <- normalizePath(args[1])
reps <- if (length(args) > 1) as.integer(args[2]) else 20
cores <- if (length(args) > 2) as.integer(args[3]) else 1
library(evomod2)

N <- 500; TRAITS <- 4; LOCI <- 50
MU <- 0.001; MU_B <- 0.0005; SIGMA <- 0.2; V_E <- 0.8
BURN_IN <- 200; STAB <- 1000; DIRECT <- 1000
DELTA <- c(0.01, 0.01, 0, 0)
OMEGA <- 10 * kronecker(diag(2), matrix(0.5, 2, 2)) + diag(5, TRAITS)  # two modules

# ---- C runs -----------------------------------------------------------------
work <- tempfile("evomod_c")
dir.create(work)
exe <- file.path(work, "evomod.exe")
cflags <- system2("gsl-config", "--cflags", stdout = TRUE)
libs <- system2("gsl-config", "--libs", stdout = TRUE)
src <- file.path(c_dir, c("main.c", "population/evolution.c", "population/moments.c", "population/io.c"))
if (system(paste("cc -O2", cflags, "-o", exe, paste(src, collapse = " "), libs)) != 0)
  stop("could not build the C simulation")

# system2 runs in the current directory; run each replicate from its own folder
run_c <- function(rep) {
  dir <- file.path(work, paste0("rep", rep))
  dir.create(file.path(dir, "input"), recursive = TRUE)
  dir.create(file.path(dir, "output"))
  write.table(OMEGA, file.path(dir, "input/omega.csv"), row.names = FALSE, col.names = FALSE)
  writeLines(paste(rep(0, TRAITS), collapse = " "), file.path(dir, "input/theta.csv"))
  writeLines(c(N, TRAITS, LOCI, MU, MU_B, SIGMA, V_E, BURN_IN, STAB, DIRECT,
               paste(DELTA, collapse = " ")), file.path(dir, "input/par.txt"))
  writeLines(c("run", "par.txt", "0", "1"), file.path(dir, "answers"))
  status <- system(sprintf("cd %s && %s < answers > log 2>&1", shQuote(dir), shQuote(exe)))
  if (status != 0) stop("C run ", rep, " failed")
  # pop.pop: per individual, y (2 x loci values) then B (traits x alleles, row major)
  raw <- scan(file.path(dir, "output/run/pop.pop"), quiet = TRUE)
  per_ind <- 2 * LOCI * (TRAITS + 1)
  x <- vapply(seq_len(N), function(i) {
    v <- raw[(i - 1) * per_ind + seq_len(per_ind)]
    y <- v[seq_len(2 * LOCI)]
    B <- matrix(v[-seq_len(2 * LOCI)], TRAITS, byrow = TRUE)
    drop(B %*% y)
  }, numeric(TRAITS))
  summarise(x)
}

summarise <- function(x) {
  G <- cov(t(x))
  c(G11 = G[1, 1], G33 = G[3, 3], G12 = G[1, 2], G34 = G[3, 4], G13 = G[1, 3],
    x1 = mean(x[1, ]), x3 = mean(x[3, ]))
}

# ---- evomod2 runs -------------------------------------------------------------
model <- evomod_model(binary_b(LOCI, TRAITS, init_b = 0.25, init_y = "mutations"),
                      binary_b_mutation(MU, SIGMA, MU_B, b_mutation = "flip"),
                      env_sd = V_E, init_sd = SIGMA)
run_r <- function(rep) {
  pop <- evolve(model, n = N, record_every = 1000, seed = rep,
                regimes = list(drift(BURN_IN), stabilizing(STAB, omega = OMEGA),
                               directional(DIRECT, delta_theta = DELTA, omega = OMEGA)))
  summarise(pop$x)
}

c_res <- do.call(rbind, parallel::mclapply(seq_len(reps), run_c, mc.cores = cores))
r_res <- do.call(rbind, parallel::mclapply(seq_len(reps), run_r, mc.cores = cores))

summary <- do.call(rbind, lapply(colnames(c_res), function(s) {
  se <- sqrt(var(r_res[, s]) / reps + var(c_res[, s]) / reps)
  data.frame(statistic = s, evomod2 = mean(r_res[, s]), C = mean(c_res[, s]),
             z = (mean(r_res[, s]) - mean(c_res[, s])) / se)
}))
print(summary, digits = 3, row.names = FALSE)
cat("\n|z| > 3 for:", if (any(abs(summary$z) > 3))
  paste(summary$statistic[abs(summary$z) > 3], collapse = ", ") else "none", "\n")
