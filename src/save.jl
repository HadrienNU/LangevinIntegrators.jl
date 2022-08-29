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
end

function TrajectorySave(n_save_iters::Int, file_pattern::Union{String,Nothing}, id_traj::Int, n_step::Int, state::AbstractState; kwargs...)
    n_save = fld(n_step , n_save_iters) # max number of step to save
    if isnothing(file_pattern)
        file=nothing
        buffer_size = n_save  # If no file, then set buffer size to save full traj
    else
        filename = replace(file_pattern,"*"=>id_traj)
        file = open(filename, "w")
        buffer_size = get(kwargs,:buffer_size, n_save)
    end

    traj = Array{typeof(state.x[1])}(undef, (buffer_size , size(state.x)[1]))
    time = Array{Float64}(1:buffer_size)
    return TrajectorySave(n_save_iters, time, traj, 0, buffer_size,file)
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

# To have output as x or (x,v) or (x,v,h) either define different struct with different array or implement just one big array

# mutable struct TrajectorySaveInertial <: AbstractSave
#     n_save_iters::Int
#     time::Array
#     xt::Array
#     vt::Array
#     save_index::Int
#     # to_save::Array{Range} # Variable that say what to save, by defaut only save position
#
#     buffer_size::Int # Par défaut on va dire que toute la traj est mise dans le fichier en une seule fois
#     file::Union{IOStream,Nothing} # When nothing, we don't save to file
# end
#
#
#
# mutable struct TrajectorySaveHidden <: AbstractSave
#     n_save_iters::Int
#     time::Array
#     xt::Array
#     vt::Array
#     ht::Array
#     save_index::Int
#     # to_save::Array{Range} # Variable that say what to save, by defaut only save position
#
#     buffer_size::Int # Par défaut on va dire que toute la traj est mise dans le fichier en une seule fois
#     file::Union{IOStream,Nothing} # When nothing, we don't save to file
# end

# function get_state_size(save::TrajectorySaveHidden,state)
#     return size(state.x)[1]+size(state.v)[1]+size(state.h)[1]
# end
#
# function get_state(save::TrajectorySaveHidden,state)
#     return hcat(deepcopy(state.x),deepcopy(state.v),deepcopy(state.h))
# end


#In case, there is no traj_save

function save_state(save::Nothing, t::Float64, state::AbstractState; kwargs...)
end

function flush_traj(save::Nothing) # At the end of the traj write to file the rest of the data

end
