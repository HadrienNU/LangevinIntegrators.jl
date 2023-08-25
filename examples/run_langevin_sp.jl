#=
run_kernel.jl

Run a simple trajectory for a kernel
=#

using LangevinIntegrators
using Random

Random.seed!(7)

force = ForceFromPotential("Harmonic")
params = TrajsParams(
    n_steps = 10^5,
    n_trajs = 1,
)
Δt = 1e-3
β = 1.0

integrator = VECSp(force, β, x->1.0, 1.0, Δt, 1) 
init_conds = initialize_initcond(integrator, Dict())
init_state = InitState(integrator, init_conds)

save_trajs = TrajectorySave(
    params.n_save_iters,
    params.save_filename_pattern,
    5,
    params.n_steps,
    init_state,
    to_save = [:x, :v_mid],
)
#save_trajs = nothing
Random.seed!(7)
run_trajectory!(init_state, integrator, save_trajs; params = params)
