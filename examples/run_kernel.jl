#=
run_kernel.jl

Run a simple trajectory for a kernel
=#

using LangevinIntegrators
using Random

Random.seed!(7)

force=ForceFromPotential("Harmonic")
params=TrajsParams(n_steps = 10^5, n_trajs = 1, n_save_iters = 50, save_filename_pattern = "*-traj.dat")
Δt = 1e-3
tps = LinRange(0,div(params.n_steps,2)*Δt,div(params.n_steps,2)+1)
kernel= exp.(-20.0*LinRange(0,500*1e-3, 500))
β = 1.0

integrator=BBK_Kernel(force, β , kernel[1:1], 1.0, Δt, 1)  #BBK_Kernel,GJF_Kernel
# integrator=BBK(force, β , kernel[1], 1.0, Δt, 1)  #BBK_Kernel,GJF_Kernel,EM_Kernel
init_conds = initialize_initcond(integrator, Dict())
init_state = InitState(integrator, init_conds)

save_trajs = TrajectorySave(params.n_save_iters, params.save_filename_pattern, 1, params.n_steps, init_state)
#save_trajs = nothing
Random.seed!(7)
run_trajectory!(init_state, integrator, save_trajs; params = params)
