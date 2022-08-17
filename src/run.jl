#=
run.jl

Launch the trajectories
=#

function init_trajectory(integrator::S; params = LangevinParams()) where {S<:AbstractIntegrator}

		state = InitState!([0.0], integrator)
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
        end
#        for observer in # On itere sur tous les Observer, qui sont soit des analyse statistique, soit des dump de la traj (en memoire ou en fichier),
#            AnalyzeTrajectory!(observer,state) # Compute observables and dump data if required
#		end
#        CheckStoppingCriterion()  # For committor and FPT computations
    end
	return traj
end

# """
# Pour lancer plein de traj en parallèle
# """
# function run_trajectories(args)
# 	body
# end
