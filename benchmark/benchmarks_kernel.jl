using BenchmarkTools
using LangevinIntegrators


#Defining somes benchmark

benchmarkdir = abspath(dirname(@__FILE__))


init_conds_args =
    Dict("position" => Dict("type" => "Cste"), "velocity" => Dict("type" => "Cste"))
force = ForceFromPotential("Harmonic")
params = TrajsParams(; n_steps = 10^5)
β = 1.0
Δt = 1e-3

integrator_em = EM(force, β, 0.001)
init_conds_em = initialize_initcond(integrator_em, init_conds_args)

integrator_bbk = BBK(force, β, 1.0, 1.0, Δt)

const SUITE = BenchmarkGroup()


SUITE["integrator_kernel"] = BenchmarkGroup(["run"])
#Benchmark for Kernel integrator
kernel = exp.(-20.0 * LinRange(0, 500 * Δt, 500))

integrator_kernelbbk = BBK_Kernel(force, β, kernel, 1.0, Δt, 1)
init_conds = initialize_initcond(integrator_kernelbbk, init_conds_args)
state_kernelbbk = InitState(integrator_kernelbbk, init_conds)
SUITE["integrator_kernel"]["run_BBK"] = @benchmarkable run_trajectory!(
    $state_kernelbbk,
    $integrator_kernelbbk;
    params = $params,
)

integrator_kernelgjf = GJF_Kernel(force, β, kernel, 1.0, Δt, 1)
init_conds = initialize_initcond(integrator_kernelgjf, init_conds_args)
state_kernelgjf = InitState(integrator_kernelgjf, init_conds)
SUITE["integrator_kernel"]["run_GJF"] = @benchmarkable run_trajectory!(
    $state_kernelgjf,
    $integrator_kernelgjf;
    params = $params,
)
