# Shared settings; keep in sync with cross_check.R
const N, LOCI, P = 500, 50, 2
const MU, MUT_SD, ENV_SD, INIT_SD = 0.005, 0.05, 1.0, 0.05
const OMEGA, STAB, DIRECT = 20.0, 1000, 1000
const DELTA = [0.005, 0.0]
self_eps(e) = (E = zeros(P, P, P); for k in 1:P; E[k, k, k] = e; end; E)
const MODELS = ["additive" => nothing, "multilinear" => self_eps(0.05)]
