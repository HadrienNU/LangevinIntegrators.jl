#=
run.jl

Launch the trajectories
=#

function init_trajectory(integrator::S; params = LangevinParams()) where {S<:AbstractIntegrator}
		state = InitState(integrator)
		#Ensuite on peut randomiser les vitesses
		return state
end


function run_trajectory!(state::IS, integrator::S; params = LangevinParams()) where {IS<:AbstractState,S<:AbstractIntegrator}
	traj = Array{typeof(state.x[1])}(undef,(params.n_save,size(state.x)[1])) #[similar(State.x) for i = 1:params.n_save]
	save_index=1
	for n = 1:params.n_iters
        	UpdateState!(state, integrator)
			if (mod(n, params.n_save_iters) == 0)
	            @. traj[save_index,:] = deepcopy(state.x)
	            save_index += 1
	        # end
	#        for observer in # On itere sur tous les Observer, qui sont soit des analyse statistique, soit des dump de la traj (en memoire ou en fichier),
	#            AnalyzeTrajectory!(observer,state) # Compute observables and dump data if required
	#		end
	#        CheckStoppingCriterion()  # For committor and FPT computations

		end
    end
	return traj
end

# """
# Pour lancer plein de traj en parallèle
# """
function run_trajectories(integrator::S; params = LangevinParams()) where {S<:AbstractIntegrator}
	state = InitState(integrator) # To get the type of the state
	init_states=Vector{typeof(state)}(undef,params.n_trajs)
	for n in 1:params.n_trajs
		init_states[n] = init_trajectory(integrator;params=params)
	end
	Threads.@threads for n in 1:params.n_trajs # If there is only only Thread that would be serial
           run_trajectory!(init_states[n], integrator; params = params)
       end
end
