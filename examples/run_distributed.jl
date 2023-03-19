#=
run_distributed.jl

An example of distributed computing of the trajectories. It should be run with julia -p [nprocs >= 2]

Res hold the return value of run_trajectory!
init_states is not changed

=#
@everywhere using LangevinIntegrators

let
    force = ForceFromPotential("Harmonic")
    integrator = EM(force, 1.0, 0.001)
    params = TrajsParams(;
        n_steps = 1e7,
        n_trajs = 20,
        n_save_iters = 2e6,
        save_filename_pattern = "trajectory_*.dat",
    )
    init_states = generate_initial_conditions(integrator; params = params)
    #Dans le cas distibué, il faut fournir des couples init_states, integrator
    run = state -> run_trajectory!(state, integrator; params = params)
    res = pmap(run, init_states)
    for n = 1:params.n_trajs
        println(res)
    end
end
