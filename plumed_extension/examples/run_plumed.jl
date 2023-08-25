#=
main.jl

Read parameters file and launch the trajectories
=#

# using Plots
using LangevinIntegrators
using LangevinIntegratorsPlumedExt

let
    force = ForceFromPotential("Flat", 6)
    integrator = Verlet(force, 1.0, 0.001, 6)
    params = TrajsParams()
    state = InitState([1.0,2.0,3.0,4.0,7.0,1.0], [1.0,2.0,3.0,4.0,7.0,1.0],integrator, )
    addPlumed!(integrator, "plumed.dat", "plumed.log")

    traj = TrajectorySave(
        params.n_save_iters,
        params.save_filename_pattern,
        1,
        params.n_steps,
        state,
    )
    run_trajectory!(state, integrator, traj; params = params)
    println(state.x)
    # plot!(traj.time, traj.xt[1][:, 1], label = "Traj")
    # xlabel!("t")
    # ylabel!("x")
end
