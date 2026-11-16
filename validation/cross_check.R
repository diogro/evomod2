# Cross-check evomod2 against the reference Julia implementation (EvoMod.jl).
#
# Runs the same model and regimes in both and compares replicate means of the
# final G matrix and mean genotypic values. Implementations use different
# random number generators, so results agree in distribution, not exactly.
#
# Usage: Rscript validation/cross_check.R <path to evomod_julia> [reps]

args <- commandArgs(trailingOnly = TRUE)
julia_pkg <- args[1]
reps <- if (length(args) > 1) as.integer(args[2]) else 20
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
library(evomod2)

# settings: keep in sync with cross_check_settings.jl
N <- 500; LOCI <- 50; P <- 2
MU <- 0.005; MUT_SD <- 0.05; ENV_SD <- 1; INIT_SD <- 0.05
OMEGA <- 20; STAB <- 1000; DIRECT <- 1000; DELTA <- c(0.005, 0)
self_eps <- function(e) { E <- array(0, c(P, P, P)); for (k in 1:P) E[k, k, k] <- e; E }
architectures <- list(additive = additive(LOCI, P), multilinear = multilinear(LOCI, P, self_eps(0.05)))

r_rows <- do.call(rbind, lapply(names(architectures), function(name) {
  model <- evomod_model(architectures[[name]], mutation(MU, MUT_SD), env_sd = ENV_SD, init_sd = INIT_SD)
  do.call(rbind, lapply(seq_len(reps), function(rep) {
    pop <- evolve(model, stabilizing(STAB, omega = OMEGA), n = N, record_every = STAB, seed = rep)
    G_stab <- cov(t(pop$x))
    pop <- evolve(pop, directional(DIRECT, delta_theta = DELTA, omega = OMEGA))
    G <- cov(t(pop$x))
    data.frame(model = name, rep = rep, G11_stab = G_stab[1, 1], G11 = G[1, 1], G22 = G[2, 2],
               G12 = G[1, 2], x1 = mean(pop$x[1, ]), x2 = mean(pop$x[2, ]))
  }))
}))

jl_out <- tempfile(fileext = ".csv")
status <- system2("julia", c(paste0("--project=", julia_pkg), "--startup-file=no",
                             file.path(here, "cross_check.jl"), reps, jl_out))
if (status != 0) stop("Julia run failed")
jl_rows <- read.csv(jl_out)

stats <- c("G11_stab", "G11", "G22", "G12", "x1", "x2")
summary <- do.call(rbind, lapply(names(architectures), function(name) {
  do.call(rbind, lapply(stats, function(s) {
    r <- r_rows[r_rows$model == name, s]; j <- jl_rows[jl_rows$model == name, s]
    se <- sqrt(var(r) / length(r) + var(j) / length(j))
    data.frame(model = name, statistic = s, R = mean(r), Julia = mean(j),
               z = (mean(r) - mean(j)) / se)
  }))
}))
print(summary, digits = 3, row.names = FALSE)
cat("\n|z| > 3 for:", if (any(abs(summary$z) > 3)) paste(summary$statistic[abs(summary$z) > 3], collapse = ", ") else "none", "\n")
