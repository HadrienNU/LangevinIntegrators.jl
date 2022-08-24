#=
run_distributed.jl

An example of distributed computing of the trajectories

In this case, we should implement a dump struct that hold the id of the initial state if we want to save to files

Res hold the return value of run_trajectory!
init_states is not changed

=#
@everywhere using LangevinIntegrators

let
    force=ForceFromPotential("Harmonic")
    integrator=EM(force,1.0,0.001)
    params=TrajsParams(;n_steps=1e7,n_trajs=20,n_save_iters=2e6)
    init_states=generate_initial_conditions(integrator; params = params)
    #Dans le cas distibué, il faut fournir des couples init_states, integrator
    run=state->run_trajectory!(state, integrator; params = params)
    res=pmap(run, init_states);
    for n in 1:params.n_trajs
        println(res[n])
    end
end
