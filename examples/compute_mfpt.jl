#=
run_kernel.jl

Run a simple trajectory for a kernel
=#

using LangevinIntegrators

force = ForceFromPotential("DoubleWell",potential_kwargs=Dict(:a=>2.0))
addFix!(force,Stopper([0.5]))
params = TrajsParams(
    n_steps = 10^5,
    n_trajs = 1000,
)
Δt = 1e-3
gamma =1.0
β = 1.0

init_conds = Dict(
    "position" => Dict("type" => "cste", "val" => -1.0),
    "velocity" => Dict("type" => "cste", "val" => 0.0),
    # "velocity" => Dict("type" => "gaussian", "std" => sqrt(1 / β)),
)

integrator = VEC(force, β,gamma, 1.0, Δt, 1)
fpt,reached =run_fpt(integrator; params = params, init_conds_args = init_conds)


println("MFPT estimate: ",sum(fpt)/sum(reached)) #," +/- ", )
