#=
main.jl

Read parameters file and launch the trajectories
=#

using Plots
using LangevinIntegrators
using LangevinIntegratorsPlumedExt

let
    force=ForceFromPotential("Flat")
    integrator=EM(force,1.0,0.001)
    params=TrajsParams()
	state=InitState(integrator, params)
    addPlumed!(integrator,"plumed.dat", "plumed.log")

    traj = TrajectorySave(params.n_save_iters, params.save_filename_pattern, 1, params.n_steps, state)
    run_trajectory!(state, integrator, traj; params =params)
    println(state.x)
    plot!(traj.time, traj.xt[1][:,1],label="Traj")
    xlabel!("t")
    ylabel!("x")
end
