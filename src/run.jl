#=
run.jl

Launch the trajectories
=#


function generate_initial_conditions(integrator::S; params = TrajsParams(), init_conds_args = Dict()) where {S<:AbstractIntegrator}
    init_conds = initialize_initcond(integrator, init_conds_args)
    state = InitState(integrator, init_conds; id = 1) # To get the type of the state
    init_states = Vector{typeof(state)}(undef, params.n_trajs)
    init_states[1] = state
    for n = 2:params.n_trajs
        init_states[n] = InitState(integrator, init_conds; id = n)
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
function run_trajectory!(state::IS, integrator::S, save_traj::Union{AbstractSave,Nothing} = nothing; params = TrajsParams(), kwargs...) where {IS<:AbstractState,S<:AbstractIntegrator}
    n = 0
    for fix in integrator.force.fixes
        init_fix(fix; kwargs...)
    end
    #Update initially force ??
    nostop = 0
    # nostop = forceUpdate!(integrator.force, state.f, state.x; step = n) # Except that depending of the integrator state.f does not exist
    while n < params.n_steps && nostop == 0
        n += 1
        nostop = UpdateState!(state, integrator; step = n)
        if !isnothing(save_traj) && (mod(n, save_traj.n_save_iters) == 0)
            save_state(save_traj, n * integrator.Δt, state; kwargs...) # Compute observables and dump data if required
        end
    end
    # End of the fix to close file if needed
    for fix in integrator.force.fixes
        close_fix(fix)
    end

    flush_traj(save_traj) # Save the remaining of the traj to file if wanted

    return n * integrator.Δt, nostop # Ca va permettre les calculs de temps de transitions facilement
end

"""
For launching a bunch of trajectories. For distributed computing see examples
"""
function run_trajectories(integrator::S; params = TrajsParams(), init_conds_args = Dict(), kwargs...) where {S<:AbstractIntegrator}
    init_states = generate_initial_conditions(integrator; params = params, init_conds_args = init_conds_args)
    #To make the system threadsafe, we have to copy the integrator into nthreads copy, and provide one copy per thread
    if Threads.nthreads() > 1
        integrators_set = [deepcopy(integrator) for n = 1:Threads.nthreads()]
    else
        integrators_set = [integrator]
    end
    #Il faudrait aussi faire un truc pour créer un dossier par thread pour écrire les fichiers pour que ça ne se marche pas dessus (notamment pour plumed)
    save_trajs = Array{AbstractSave}(undef,params.n_trajs)
    Threads.@threads for n = 1:params.n_trajs # If there is only only Thread that would be serial
        save_trajs[n] = TrajectorySave(params.n_save_iters, params.save_filename_pattern, n, params.n_steps, init_states[n]; kwargs...)
        run_trajectory!(init_states[n], integrators_set[Threads.threadid()], save_trajs[n]; params = params, id_traj=n, kwargs...)
    end
    return save_trajs # This is a set of trajectories
end


function run_fpt(integrator::S; params = TrajsParams(), init_conds_args = Dict(), kwargs...) where {S<:AbstractIntegrator}
    init_states = generate_initial_conditions(integrator; params = params, init_conds_args = init_conds_args)
    #To make the system threadsafe, we have to copy the integrator into nthreads copy, and provide one copy per thread
    if Threads.nthreads() > 1
        integrators_set = [deepcopy(integrator) for n = 1:Threads.nthreads()]
    else
        integrators_set = [integrator]
    end

    reached = Array{Bool}(undef,params.n_trajs)
    fpt=Array{Float64}(undef,params.n_trajs)

    Threads.@threads for n = 1:params.n_trajs # If there is only only Thread that would be serial
        finaltime, stoppingCriterion = run_trajectory!(init_states[n], integrators_set[Threads.threadid()], nothing; params = params, id_traj=n, kwargs...)
        reached[n] = (stoppingCriterion != 0)
        fpt[n] = finaltime
    end
    return fpt, reached # This is a set of trajectories
end
