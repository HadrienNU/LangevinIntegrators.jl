#=
save.jl

The code to save the generated trajectories
=#

# The code to save the trajectory, either in an array, or in a file
# For now only save the position
# To save more, maybe juste create one struct per type of state and save everything that we can
# If more specific saving is needed for performance reasons, then we can just define a custom save taylord for this
# TODO: write also a xyz file format type save, that could be used if file_pattern end with xyz

abstract type AbstractSave end

mutable struct TrajectorySave <: AbstractSave
    n_save_iters::Int
    time::Array
    xt::Array
    save_index::Int
    # to_save::Array{Range} # Variable that say what to save, by defaut only save position

    buffer_size::Int # Par défaut on va dire que toute la traj est mise dans le fichier en une seule fois
    file::Union{IOStream,Nothing} # When nothing, we don't save to file
    function TrajectorySave(n_save_iters::Int, filename::Union{String,Nothing}, buffer_size::Int, state::AbstractState)
        file = isnothing(filename) ? nothing : open(filename, "w")
        time = Array{Float64}(1:buffer_size)
        xt = Array{typeof(state.x[1])}(undef, (buffer_size , size(state.x)[1]))
        new(n_save_iters, time, xt, 0, buffer_size, file)
    end
end

function TrajectorySave(n_save_iters::Int, file_pattern::Union{String,Nothing}, id_traj::Int, n_step::Int, state::AbstractState; kwargs...)
    n_save = fld(n_step , n_save_iters) # max number of step to save
    if isnothing(file_pattern)
        filename=nothing
        buffer_size = n_save  # If no file, then set buffer size to save full traj
    else
        filename = replace(file_pattern,"*"=>id_traj)
        buffer_size = get(kwargs,:buffer_size, n_save)
    end
    # Then we can select the save struct based on kwargs
    save_velocity = typeof(state) <: AbstractInertialState && get(kwargs,:save_velocity, false)
    save_hidden =  (typeof(state)<: AbstractOverdampedMemoryHiddenState || typeof(state)<: AbstractMemoryHiddenState) && get(kwargs,:save_hidden, false)

    if save_velocity && save_hidden
        return TrajectorySaveHidden(n_save_iters, filename, buffer_size, state)
    elseif save_velocity && !save_hidden
        return TrajectorySaveInertial(n_save_iters, filename, buffer_size, state)
    elseif save_hidden
        return TrajectorySaveOnlyHidden(n_save_iters, filename, buffer_size, state)
    end
    return TrajectorySave(n_save_iters, filename, buffer_size, state)
end


function save_state(save::TrajectorySave, t::Float64, state::AbstractState; kwargs...)
    save.save_index += 1
    save.xt[save.save_index, :] .= deepcopy(state.x)
    save.time[save.save_index] = t
    # On peut faire quelque pour enregistrer dans un fichier tous les nb pas de temps, on remet alors save_index à 1
    if  !isnothing(save.file) && save.save_index >= save.buffer_size
        writedlm(save.file, hcat(save.time[1:save.save_index], save.xt[1:save.save_index,:]))
        save.save_index = 0
    end

end

function flush_traj(save::TrajectorySave) # At the end of the traj write to file the rest of the data
    if !isnothing(save.file)
        writedlm(save.file, hcat(save.time[1:save.save_index], save.xt[1:save.save_index,:]))
        close(save.file)
    end
end


mutable struct TrajectorySaveInertial <: AbstractSave
    n_save_iters::Int
    time::Array
    xt::Array
    vt::Array
    save_index::Int

    buffer_size::Int # Par défaut on va dire que toute la traj est mise dans le fichier en une seule fois
    file::Union{IOStream,Nothing} # When nothing, we don't save to file

    function TrajectorySaveInertial(n_save_iters::Int, filename::Union{String,Nothing}, buffer_size::Int, state::AbstractInertialState)
        file = isnothing(filename) ? nothing : open(filename, "w")
        time = Array{Float64}(1:buffer_size)
        xt = Array{typeof(state.x[1])}(undef, (buffer_size , size(state.x)[1]))
        vt = Array{typeof(state.v[1])}(undef, (buffer_size , size(state.v)[1]))
        new(n_save_iters, time, xt, vt, 0, buffer_size, file)
    end
end

function save_state(save::TrajectorySaveInertial, t::Float64, state::AbstractState; kwargs...)
    save.save_index += 1
    save.xt[save.save_index, :] .= deepcopy(state.x)
    save.vt[save.save_index, :] .= deepcopy(state.v)
    save.time[save.save_index] = t
    # On peut faire quelque pour enregistrer dans un fichier tous les nb pas de temps, on remet alors save_index à 1
    if  !isnothing(save.file) && save.save_index >= save.buffer_size
        writedlm(save.file, hcat(save.time[1:save.save_index], save.xt[1:save.save_index,:], save.vt[1:save.save_index,:]))
        save.save_index = 0
    end

end

function flush_traj(save::TrajectorySaveInertial) # At the end of the traj write to file the rest of the data
    if !isnothing(save.file)
        writedlm(save.file, hcat(save.time[1:save.save_index], save.xt[1:save.save_index,:], save.vt[1:save.save_index,:]))
        close(save.file)
    end
end


mutable struct TrajectorySaveOnlyHidden <: AbstractSave
    n_save_iters::Int
    time::Array
    xt::Array
    ht::Array
    save_index::Int

    buffer_size::Int # Par défaut on va dire que toute la traj est mise dans le fichier en une seule fois
    file::Union{IOStream,Nothing} # When nothing, we don't save to file

    function TrajectorySaveOnlyHidden(n_save_iters::Int, filename::Union{String,Nothing}, buffer_size::Int, state::Union{AbstractOverdampedMemoryHiddenState,AbstractMemoryHiddenState})
        file = isnothing(filename) ? nothing : open(filename, "w")
        time = Array{Float64}(1:buffer_size)
        xt = Array{typeof(state.x[1])}(undef, (buffer_size , size(state.x)[1]))
        ht = Array{typeof(state.h[1])}(undef, (buffer_size , size(state.h)[1]))
        new(n_save_iters, time, xt, ht, 0, buffer_size, file)
    end
end


function save_state(save::TrajectorySaveOnlyHidden, t::Float64, state::AbstractState; kwargs...)
    save.save_index += 1
    save.xt[save.save_index, :] .= deepcopy(state.x)
    save.ht[save.save_index, :] .= deepcopy(state.h)
    save.time[save.save_index] = t
    # On peut faire quelque pour enregistrer dans un fichier tous les nb pas de temps, on remet alors save_index à 1
    if  !isnothing(save.file) && save.save_index >= save.buffer_size
        writedlm(save.file, hcat(save.time[1:save.save_index], save.xt[1:save.save_index,:], save.ht[1:save.save_index,:]))
        save.save_index = 0
    end

end

function flush_traj(save::TrajectorySaveOnlyHidden) # At the end of the traj write to file the rest of the data
    if !isnothing(save.file)
        writedlm(save.file, hcat(save.time[1:save.save_index], save.xt[1:save.save_index,:], save.ht[1:save.save_index,:]))
        close(save.file)
    end
end

mutable struct TrajectorySaveHidden <: AbstractSave
    n_save_iters::Int
    time::Array
    xt::Array
    vt::Array
    ht::Array
    save_index::Int

    buffer_size::Int # Par défaut on va dire que toute la traj est mise dans le fichier en une seule fois
    file::Union{IOStream,Nothing} # When nothing, we don't save to file

    function TrajectorySaveHidden(n_save_iters::Int, filename::Union{String,Nothing}, buffer_size::Int, state::AbstractMemoryHiddenState)
        file = isnothing(filename) ? nothing : open(filename, "w")
        time = Array{Float64}(1:buffer_size)
        xt = Array{typeof(state.x[1])}(undef, (buffer_size , size(state.x)[1]))
        vt = Array{typeof(state.v[1])}(undef, (buffer_size , size(state.v)[1]))
        ht = Array{typeof(state.h[1])}(undef, (buffer_size , size(state.h)[1]))
        new(n_save_iters, time, xt, vt, ht, 0, buffer_size, file)
    end
end

function save_state(save::TrajectorySaveHidden, t::Float64, state::AbstractState; kwargs...)
    save.save_index += 1
    save.xt[save.save_index, :] .= deepcopy(state.x)
    save.vt[save.save_index, :] .= deepcopy(state.v)
    save.ht[save.save_index, :] .= deepcopy(state.h)
    save.time[save.save_index] = t
    # On peut faire quelque pour enregistrer dans un fichier tous les nb pas de temps, on remet alors save_index à 1
    if  !isnothing(save.file) && save.save_index >= save.buffer_size
        writedlm(save.file, hcat(save.time[1:save.save_index], save.xt[1:save.save_index,:], save.vt[1:save.save_index,:], save.ht[1:save.save_index,:]))
        save.save_index = 0
    end

end

function flush_traj(save::TrajectorySaveHidden) # At the end of the traj write to file the rest of the data
    if !isnothing(save.file)
        writedlm(save.file, hcat(save.time[1:save.save_index], save.xt[1:save.save_index,:], save.vt[1:save.save_index,:], save.ht[1:save.save_index,:]))
        close(save.file)
    end
end




mutable struct TransitionObserver <: AbstractSave
    n_save_iters::Int
    x0::Float64
    x1::Float64
    timesAB::Array
    timesBA::Array

    curr_ind::Int
    last_ind::Int
    former_ind::Int

    last_crossing_time::Float64
    curr_transition_times::Float64

    function TransitionObserver(n_save_iters::Int, x0::Float64,x1::Float64, state::AbstractState; kwargs...)
        timesAB = Array{Float64,1}()
        timesBA = Array{Float64,1}()

        ind = state.x[1] <= x0 ? -1 : (state.x[1] >= x1 ? 1 : 0)

        new(n_save_iters, x0, x1, timesAB, timesBA, ind, ind, 0, 0.0, 0.0)
    end
end


function save_state(save::TransitionObserver, t::Float64, state::AbstractState; kwargs...)

    #D'abord on récupére l'index de l'état actuel
    save.curr_ind = state.x[1] <= save.x0 ? -1 : (state.x[1] >= save.x1 ? 1 : 0)
    if save.curr_ind != save.last_ind
        save.curr_transition_times += t-save.last_crossing_time
        save.last_crossing_time = t
        if save.curr_ind != save.former_ind # If we are doing -1 0 1 or 1 0 -1
            if save.curr_ind == -1
                append!(save.timesBA,save.curr_transition_times)
            elseif save.curr_ind == 1
                append!(save.timesAB,save.curr_transition_times)
            end
            save.curr_transition_times = 0.0
        end
        save.former_ind = save.last_ind
        save.last_ind = save.curr_ind
    end
end

function flush_traj(save::TransitionObserver) # At the end of the traj write to file the rest of the data

end

#In case, there is no traj_save

function save_state(save::Nothing, t::Float64, state::AbstractState; kwargs...)
end

function flush_traj(save::Nothing) # At the end of the traj write to file the rest of the data

end
