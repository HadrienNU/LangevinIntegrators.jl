#=
run_kernel.jl

Run a simple trajectory for a kernel
=#

using LangevinIntegrators

force=ForceFromPotential("Harmonic")
params=TrajsParams(n_steps = 10^4, n_trajs = 1, n_save_iters = 50)
Δt = 1e-3
tps = LinRange(0,div(params.n_steps,2)*Δt,div(params.n_steps,2)+1)
kernel= exp.(-20.0*LinRange(0,500*1e-3, 500))
β = 1.0

integrator=GJF_Kernel(force, β , kernel, 1.0, Δt, 1)  #BBK_Kernel,GJF_Kernel,EM_Kernel
init_conds = initialize_initcond(integrator, Dict())
init_state = InitState(integrator, init_conds)

# save_trajs = TrajectorySave(params.n_save_iters, params.save_filename_pattern, 1, params.n_steps, init_states; kwargs...)
save_trajs = nothing
run_trajectory!(init_state, integrator, save_trajs; params = params)
