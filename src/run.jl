#=
run.jl

Launch the trajectories
=#

function generate_initial_conditions(integrator::S; params = LangevinParams(), init_conds_args=Dict()) where {S<:AbstractIntegrator}
    init_conds=initialize_initcond(integrator;init_conds_args...)
	state = InitState(integrator) # To get the type of the state
	init_states=Vector{typeof(state)}(undef,params.n_trajs)
	for n in 1:params.n_trajs
		init_states[n] = InitState(integrator,init_conds;id=n)
	end
    return init_states
end


"""
    run_trajectory!(state, integrator)
evolve in time one trajectory starting at state and integrated by integrator
### Fields
* state   - Initial State
* integrator    - Friction matrix
"""
function run_trajectory!(state::IS, integrator::S; params = LangevinParams(), kwargs...) where {IS<:AbstractState,S<:AbstractIntegrator}
	n = 1
	nostop=0
    # Ici ik faut faire un init_fix
    for fix in integrator.force.fixes
        init_fix(fix)
    end
	while n < params.n_steps && nostop ==0
    	nostop=UpdateState!(state, integrator)
        for observer in params.observers# On itere sur tous les Observer, qui sont soit des analyse statistique, soit des dump de la traj (en memoire ou en fichier),
			if (mod(n, observer.n_save_iters) == 0)
				run_obs(observer,n*integrator.Δt,state; kwargs) # Compute observables and dump data if required
			end
		end
    end
    # Et là un end_fix
    for fix in integrator.force.fixes
        close_fix(fix)
    end
    for observer in params.observers# On itere sur tous les Observer, qui sont soit des analyse statistique, soit des dump de la traj (en memoire ou en fichier),
        end_obs_traj(observer;end_time=n*integrator.Δt,stoppingCriterion=nostop) # Close file if needed
    end
	return n*integrator.Δt,nostop # Ca va permettre les calculs de temps de transitions facilement
end

"""
For launching a bunch of trajectories. For distributed computing see examples
"""
function run_trajectories(integrator::S; params = LangevinParams(), init_conds_args=Dict(), kwargs...) where {S<:AbstractIntegrator}
    init_states=generate_initial_conditions(integrator; params = params,init_conds_args=init_conds_args)
    #To make the system threadsafe, we have to copy the integrator into nthreads copy, and provide one copy per thread
    if Threads.nthreads() > 1
        integrators_set = [deepcopy(integrator) for n in 1:Threads.nthreads()]
    else
        integrators_set = [integrator]
    end
    #Il faudrait aussi faire un truc pour créer un dossier par thread pour écrire les fichiers pour que ça ne se marche pas dessus (notamment pour plumed)

	Threads.@threads for n in 1:params.n_trajs # If there is only only Thread that would be serial
           run_trajectory!(init_states[n], integrators_set[Threads.threadid()]; params = params,kwargs...)
       end

       for observer in params.observers# On itere sur tous les Observer, qui sont soit des analyse statistique, soit des dump de la traj (en memoire ou en fichier),
           close_obs(observer) # Conclude global analysis
       end
end
