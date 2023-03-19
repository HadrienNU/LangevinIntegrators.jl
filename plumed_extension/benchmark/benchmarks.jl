using BenchmarkTools
using LangevinIntegrators
using LangevinIntegratorsPlumedExt

#Defining somes benchmark

benchmarkdir = abspath(dirname(@__FILE__))

init_conds_args =
    Dict("position" => Dict("type" => "Cste"), "velocity" => Dict("type" => "Cste"))
params = TrajsParams(; n_steps = 10^5)
β = 1.0
Δt = 1e-3

init_conds_em = initialize_initcond(integrator_em, init_conds_args)

# Benchmark with plumed
plumedFix = plumed(
    joinpath(benchmarkdir, "../test/test_plumed.dat"),
    "plumed.log",
    1,
    integrator_em.Δt;
    temperature = 1.0,
)
force_plumed = ForceFromPotential("Harmonic")
addFix!(force_plumed, plumedFix)
integrator_plumed = EM(force_plumed, 1.0, integrator_em.Δt)
state_em = InitState(integrator_plumed, init_conds_em)
SUITE["integrator_overdamped"]["run_EM_w_plumed"] =
    @benchmarkable run_trajectory!($state_em, $integrator_plumed; params = $params)
params_x10 = TrajsParams(; n_steps = 10^6)
SUITE["integrator_overdamped"]["run_EM_w_plumed_10x"] =
    @benchmarkable run_trajectory!($state_em, $integrator_plumed; params = $params_x10)
rm("plumed.log")
