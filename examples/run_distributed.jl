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
    params=LangevinParams(;n_steps=1e7,n_trajs=20,n_save_iters=2e6)
    state = InitState(integrator) # To get the type of the state
	init_states=Vector{typeof(state)}(undef,params.n_trajs)
	for n in 1:params.n_trajs
		init_states[n] = InitState(integrator,params)
	end
    run=state->run_trajectory!(state, integrator; params = params)
    res=pmap(run, init_states);
    for n in 1:params.n_trajs
        println(res[n])
    end
end
