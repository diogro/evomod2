# Replicate runs of the reference Julia implementation (EvoMod.jl,
# https://github.com/diogro/evomod_julia) for cross_check.R.
# Usage: julia --project=<path to evomod_julia> cross_check.jl <reps> <out.csv>
using EvoMod, EvoMod.Multilinear, LinearAlgebra, Random, Statistics, DelimitedFiles

include("cross_check_settings.jl")

reps, out = parse(Int, ARGS[1]), ARGS[2]
rows = Any[["model" "rep" "G11_stab" "G11" "G22" "G12" "x1" "x2"]]
for (name, eps) in MODELS, rep in 1:reps
    Random.seed!(rep)
    par = Params(LOCI, P; mu = MU, mutation_cov = MUT_SD, v_e = ENV_SD, epsilon = eps,
                 init_sigma = INIT_SD)
    pop = run_pop(par, N, ["Stab"], [STAB]; omega = OMEGA * Matrix(1.0I, P, P), progress = false,
                  thin = STAB + 1)
    xs = reduce(vcat, (ind.x' for ind in pop.pop))
    G_stab = cov(xs)
    pop = run_pop(pop, ["Direct"], [DIRECT]; delta_theta = DELTA, omega = OMEGA * Matrix(1.0I, P, P),
                  progress = false, thin = DIRECT + 1)
    xs = reduce(vcat, (ind.x' for ind in pop.pop))
    G = cov(xs); mx = vec(mean(xs, dims = 1))
    push!(rows, [name rep G_stab[1, 1] G[1, 1] G[2, 2] G[1, 2] mx[1] mx[2]])
end
writedlm(out, reduce(vcat, rows), ',')
