# Cross-check evomod2 against its Julia port, EvoMod.Sim (in EvoMod.jl), for
# the components that only exist in both since the port: binary B with the
# Moran life cycle, the contribution matrix under a fluctuating optimum,
# pairwise multilinear epistasis with viability selection, and house of cards
# mutation on a genetic map under directional selection. Settings must match
# cross_check_sim.jl. The implementations use different random number
# generators, so results agree in distribution, not exactly.
#
# Usage: Rscript validation/cross_check_sim.R <path to evomod_julia> [reps] [cores]

args <- commandArgs(trailingOnly = TRUE)
julia_pkg <- args[1]
reps <- if (length(args) > 1) as.integer(args[2]) else 30
cores <- if (length(args) > 2) as.integer(args[3]) else 1
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
library(evomod2)

N <- 300; GENS <- 1500
omega_corr <- 10 * matrix(c(1, 0.9, 0.9, 1), 2)
configs <- list(
  binary_b_moran = list(
    make = function() evomod_model(binary_b(20, 2), binary_b_mutation(0.002, 0.2, 0.001),
                                   env_sd = 0.5, init_sd = 0.2, life_cycle = moran()),
    regimes = list(stabilizing(GENS, omega = 10), directional(500, delta_theta = c(0.004, 0), omega = 10))),
  contribution_fluctuating = list(
    make = function() evomod_model(contribution_matrix(2), contribution_mutation(0.005, 0.002, 0.3),
                                   env_sd = 0.2),
    regimes = list(fluctuating(GENS, every = 10, step = 0.2, omega = 1, structure = "correlated"))),
  pairwise_viability = list(
    make = function() evomod_model(multilinear(20, 2, pairwise_epsilon(20, 2, sd = 1)),
                                   mutation(0.001, 0.2), env_sd = 1, life_cycle = viability_selection()),
    regimes = list(stabilizing(GENS, omega = omega_corr))),
  hoc_map_directional = list(
    make = function() evomod_model(additive(20, 2), hoc_mutation(0.002, 0.3), env_sd = 0.5,
                                   recombination = genetic_map(2, 0.5)),
    regimes = list(stabilizing(GENS, omega = 10), directional(500, delta_theta = c(0.01, 0), omega = 10))))

r_rows <- do.call(rbind, lapply(names(configs), function(name) {
  do.call(rbind, parallel::mclapply(seq_len(reps), function(rep) {
    set.seed(rep)
    cfg <- configs[[name]]
    pop <- evolve(cfg$make(), cfg$regimes, n = N, record_every = 10000)
    G <- cov(t(pop$x))
    data.frame(config = name, rep = rep, G11 = G[1, 1], G22 = G[2, 2], G12 = G[1, 2],
               x1 = mean(pop$x[1, ]), x2 = mean(pop$x[2, ]))
  }, mc.cores = cores))
}))

jl_out <- tempfile(fileext = ".csv")
if (system2("julia", c(paste0("--project=", julia_pkg), "--startup-file=no",
                       file.path(here, "cross_check_sim.jl"), reps, jl_out)) != 0)
  stop("Julia run failed")
jl_rows <- read.csv(jl_out)

stats <- c("G11", "G22", "G12", "x1", "x2")
summary <- do.call(rbind, lapply(names(configs), function(name) do.call(rbind, lapply(stats, function(s) {
  r <- r_rows[r_rows$config == name, s]; j <- jl_rows[jl_rows$config == name, s]
  data.frame(config = name, statistic = s, R = mean(r), Julia = mean(j),
             z = (mean(r) - mean(j)) / sqrt(var(r) / length(r) + var(j) / length(j)))
}))))
print(summary, digits = 3, row.names = FALSE)
cat("\n|z| > 3 for:", if (any(abs(summary$z) > 3))
  paste(summary$config[abs(summary$z) > 3], summary$statistic[abs(summary$z) > 3], collapse = "; ") else "none", "\n")
