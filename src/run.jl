#=
run.jl

Launch the trajectories
=#

# The code to save the trajectory, either in an array, or in a file
# For now only save the position
# To save more, maybe juste create one struct per type of state and save everything that we can
# If more specific saving is needed for performance reasons, then we can just define a custom save taylord for this
# TODO: write also a xyz file format type save, that could be used if file_pattern end with xyz

mutable struct TrajectorySave
    n_save_iters::Int
    time::Array
    xt::Array
    save_index::Int
    # to_save::Array{Range} # Variable that say what to save, by defaut only save position

    buffer_size::Int # Par défaut on va dire que toute la traj est mise dans le fichier en une seule fois
    file::Union{IOStream,Nothing} # When nothing, we don't save to file
end

function TrajectorySave(n_save_iters, file_pattern, id_traj, n_step, state; kwargs...)
    if isnothing(file_pattern)
        file=nothing
    else
        filename = replace(file_pattern,"*"=>id_traj)
        file = open(filename, "w")
    end
    n_save = floor(Int, n_step / n_save_iters) # max number of step to save
    buffer_size = get(kwargs,:buffer_size, n_save)

    traj = Array{typeof(state.x[1])}(undef, size(state.x)[1])
    time = Array{Float64}(1:n_save)
    return TrajectorySave(n_save_iters, time, traj, 1, buffer_size,file)
end


function save_state(save::TrajectorySave, t::Float64, state::AbstractState; kwargs...)
    save.xt[save.save_index, :] .= deepcopy(state.x)
    save.time[save.save_index] = t
    save.save_index += 1
    # On peut faire quelque pour enregistrer dans un fichier tous les nb pas de temps, on remet alors save_index à 1
    if  !isnothing(file) && save.save_index > save.buffer_size
        writedlm(saver.file, vcat(saver.time, saver.xt)')
        save.save_index = 1
    end
end

function flush_traj(save::TrajectorySave) # At the end of the traj write to file the rest of the data
    if !isnothing(save.file)
        writedlm(save.file, vcat(save.time[1:save.save_index], save.xt[1:save.save_index,:])')
        close(file)
    end
end

#In case, there is no traj_save

function save_state(save::Nothing, t::Float64, state::AbstractState; kwargs...)
end

function flush_traj(save::Nothing) # At the end of the traj write to file the rest of the data

end

function generate_initial_conditions(integrator::S; params = TrajsParams(), init_conds_args = Dict()) where {S<:AbstractIntegrator}
    init_conds = initialize_initcond(integrator, init_conds_args)
    state = InitState(integrator) # To get the type of the state
    init_states = Vector{typeof(state)}(undef, params.n_trajs)
    for n = 1:params.n_trajs
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
function run_trajectory!(state::IS, integrator::S, save_traj::Union{TrajectorySave,Nothing} = nothing; params = TrajsParams(), kwargs...) where {IS<:AbstractState,S<:AbstractIntegrator}
    n = 0
    for fix in integrator.force.fixes
        init_fix(fix,kwargs)
    end
    #Update initially force ??
    nostop = 0
    # nostop = forceUpdate!(integrator.force, state.f, state.x; step = n) # Except that depending of the integrator state.f does not exist
    while n < params.n_steps && nostop == 0
        n += 1
        nostop = UpdateState!(state, integrator; step = n)
        if (mod(n, save_traj.n_save_iters) == 0)
            save_state(save_traj, n * integrator.Δt, state; kwargs) # Compute observables and dump data if required
        end
    end
    # End of the fix to close file if needed
    for fix in integrator.force.fixes
        close_fix(fix)
    end
    # A la fin on sauvegarde ou non la trajectoire dans un fichier
    flush_traj(save_traj)

    # If the following is useful only for MFPT and committor, create a specific run_fpt and run_committor function
    for observer in params.observers # We iterate over observers that analyze end of traj condition
        obs_end_traj(observer, state; end_time = n * integrator.Δt, stoppingCriterion = nostop)
    end
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
    save_trajs = Array{TrajectorySave}(undef,params.n_trajs)
    Threads.@threads for n = 1:params.n_trajs # If there is only only Thread that would be serial
        save_traj[n] = TrajectorySave(params.n_save_iters, params.file_pattern, n, params.n_steps, init_states[n]; kwargs)
        run_trajectory!(init_states[n], integrators_set[Threads.threadid()], save_traj[n]; params = params, id_traj=n, kwargs...)
    end
    return save_trajs # This is a set of trajectories
end
