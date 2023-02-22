#=
save.jl

The code to save the generated trajectories
=#

# The code to save the trajectory, either in an array, or in a file
# TODO: write also a xyz file format type save, that could be used if file_pattern end with xyz

abstract type AbstractSave end

mutable struct TrajectorySave <: AbstractSave
    n_save_iters::Int
    time::Array
    xt::Array{Array}
    save_index::Int
    to_save::Array{Symbol} # Variable that say what to save, by defaut only save position

    buffer_size::Int # Par défaut on va dire que toute la traj est mise dans le fichier en une seule fois
    file::Union{IOStream,Nothing} # When nothing, we don't save to file
    function TrajectorySave(n_save_iters::Int, filename::Union{String,Nothing}, buffer_size::Int, state::AbstractState, to_save::Array{Symbol})
        file = isnothing(filename) ? nothing : open(filename, "w")
        time = Array{Float64}(1:buffer_size)
        xt = Array{Array{typeof(state.x[1])}}(undef,length(to_save))
        for (i,s) in enumerate(to_save)
            xt[i] = Array{typeof(state.x[1])}(undef, (buffer_size , size(getfield(state,s),1)))
        end
        if  !isnothing(file) # Write symbol name in top of the file
            write(file,join(vcat(["t"],[string(s) for s in to_save],["\n"]),"\t"))
        end
        new(n_save_iters, time, xt, 0, to_save,buffer_size, file)
    end
end

#To change, add a vector of Symbol that are the list of field to save, we can just add to the vector in the order wanted the field to save

function TrajectorySave(n_save_iters::Int, file_pattern::Union{String,Nothing}, id_traj::Int, n_step::Int, state::AbstractState; to_save::Union{Nothing,Array} = nothing, kwargs...)
    n_save = fld(n_step , n_save_iters) # max number of step to save
    if isnothing(file_pattern)
        filename=nothing
        buffer_size = n_save  # If no file, then set buffer size to save full traj
    else
        filename = replace(file_pattern,"*"=>id_traj)
        buffer_size = get(kwargs,:buffer_size, n_save)
    end

    # Create array of symbol with the list of field to save
    if isnothing(to_save)
        to_save = [:x]
        if get(kwargs,:save_velocity, false)
            push!(to_save,:v)
        end
        if get(kwargs,:save_hidden, false)
            push!(to_save,:h)
        end
    end
    to_save_symbol = [Symbol(s) for s in to_save if Symbol(s) in propertynames(state)] # Conversion to symbol after checking  that the field exist in the state
    if isempty(to_save_symbol)
        println("TrajectorySave: no field are marked to be saved. Please check options")
    end
    return TrajectorySave(n_save_iters, filename, buffer_size, state, to_save_symbol)
end


function save_state(save::TrajectorySave, t::Float64, state::AbstractState; kwargs...)
    save.save_index += 1
    for (i,s) in enumerate(save.to_save)
        save.xt[i][save.save_index,: ] = deepcopy(getfield(state,s))
    end
    save.time[save.save_index] = t
    # On peut faire quelque pour enregistrer dans un fichier tous les nb pas de temps, on remet alors save_index à 1
    if  !isnothing(save.file) && save.save_index >= save.buffer_size
        writedlm(save.file, hcat(save.time[1:save.save_index],[xfield[1:save.save_index,:] for xfield in save.xt]...))
        save.save_index = 0
    end

end

function flush_traj(save::TrajectorySave) # At the end of the traj write to file the rest of the data
    if !isnothing(save.file)
        writedlm(save.file, hcat(save.time[1:save.save_index], [xfield[1:save.save_index,:] for xfield in save.xt]...))
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

    burnout_time::Float64

    function TransitionObserver(n_save_iters::Int, x0::Float64,x1::Float64, state::AbstractState; kwargs...)
        timesAB = Array{Float64,1}()
        timesBA = Array{Float64,1}()

        ind = state.x[1] <= x0 ? -1 : (state.x[1] >= x1 ? 1 : 0)

        new(n_save_iters, x0, x1, timesAB, timesBA, ind, ind, 0, get(kwargs,:burnout,0.0), 0.0, get(kwargs,:burnout,0.0))
    end
end


function save_state(save::TransitionObserver, t::Float64, state::AbstractState; kwargs...)

    #D'abord on récupére l'index de l'état actuel
    save.curr_ind = state.x[1] <= save.x0 ? -1 : (state.x[1] >= save.x1 ? 1 : 0)
    if save.curr_ind != save.last_ind && t >= save.burnout_time
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
    elseif save.curr_ind != save.last_ind && t < save.burnout_time
        save.former_ind = save.last_ind
        save.last_ind = save.curr_ind
    end
end

function flush_traj(save::TransitionObserver)

end

#In case, there is no traj_save

function save_state(save::Nothing, t::Float64, state::AbstractState; kwargs...)
end

function flush_traj(save::Nothing) # At the end of the traj write to file the rest of the data

end
