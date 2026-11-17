# Replicate runs of EvoMod.Sim (the Julia port of evomod2) for cross_check_sim.R.
# Usage: julia --project=<path to evomod_julia> cross_check_sim.jl <reps> <out.csv>
using EvoMod.Sim, LinearAlgebra, Random, Statistics, DelimitedFiles

reps, out = parse(Int, ARGS[1]), ARGS[2]
N, GENS = 300, 1500
omega_corr = 10 .* [1 0.9; 0.9 1]
configs = Dict(
  "binary_b_moran" => (() -> evomod_model(binary_b(20, 2), binary_b_mutation(0.002, 0.2, 0.001);
                                          env_sd = 0.5, init_sd = 0.2, life_cycle = moran()),
                       [stabilizing(GENS; omega = 10), directional(500; delta_theta = [0.004, 0.0], omega = 10)]),
  "contribution_fluctuating" => (() -> evomod_model(contribution_matrix(2), contribution_mutation(0.005, 0.002, 0.3);
                                                    env_sd = 0.2),
                                 [fluctuating(GENS; every = 10, step = 0.2, omega = 1, structure = :correlated)]),
  "pairwise_viability" => (() -> evomod_model(multilinear(20, 2, pairwise_epsilon(20, 2; sd = 1)),
                                              mutation(0.001, 0.2); env_sd = 1, life_cycle = viability_selection()),
                           [stabilizing(GENS; omega = omega_corr)]),
  "hoc_map_directional" => (() -> evomod_model(additive(20, 2), hoc_mutation(0.002, 0.3); env_sd = 0.5,
                                               recombination = genetic_map(2, 0.5)),
                            [stabilizing(GENS; omega = 10), directional(500; delta_theta = [0.01, 0.0], omega = 10)]))
rows = Any[["config" "rep" "G11" "G22" "G12" "x1" "x2"]]
for name in sort(collect(keys(configs))), rep in 1:reps
    make, regimes = configs[name]
    Random.seed!(rep)
    pop = evolve(make(), regimes; n = N, record_every = 10000)
    G = cov(pop.x'); mx = vec(mean(pop.x, dims = 2))
    push!(rows, [name rep G[1, 1] G[2, 2] G[1, 2] mx[1] mx[2]])
end
writedlm(out, reduce(vcat, rows), ',')
