#=
run.jl

Launch the trajectories
=#

"""
    run_trajectory!(state, integrator)
evolve in time one trajectory starting at state and integrated by integrator
### Fields
* state   - In place gradient of the potential
* integrator    - Friction matrix
"""
function run_trajectory!(state::IS, integrator::S; params = LangevinParams(), kwargs...) where {IS<:AbstractState,S<:AbstractIntegrator}
	n = 1
	nostop=0
	while n < params.n_steps && nostop ==0
    	nostop=UpdateState!(state, integrator)
		for observer in params.observers# On itere sur tous les Observer, qui sont soit des analyse statistique, soit des dump de la traj (en memoire ou en fichier),
			if (mod(n, observer.n_save_iters) == 0)
				run_obs(observer,n*integrator.Δt,state; kwargs) # Compute observables and dump data if required
			end
		end
    end
	return n*integrator.Δt,nostop # Ca va permettre les calculs de temps de transitions facilement
end

"""
For launching a bunch of trajectories. For distributed computing see examples
"""
function run_trajectories(integrator::S; params = LangevinParams(); init_conds_args=Dict()) where {S<:AbstractIntegrator}
	init_conds=initialize_initcond(integrator;init_conds_args...) 
	state = InitState(integrator) # To get the type of the state
	init_states=Vector{typeof(state)}(undef,params.n_trajs)
	for n in 1:params.n_trajs
		init_states[n] = InitState(integrator,init_conds;id=n)
	end
	Threads.@threads for n in 1:params.n_trajs # If there is only only Thread that would be serial
           run_trajectory!(init_states[n], integrator; params = params)
       end
end
