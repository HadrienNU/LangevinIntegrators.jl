#=
Les functions pour le dump
=#

abstract type AbstractObserver end

abstract type AbstractDump <: AbstractObserver end
abstract type AbstractStatisticalObs <: AbstractObserver end

"""
Function to initialize the observers
"""
function setup_observers(args)
    # Ca prend une liste de dict et ça génère une observable par element de la liste
    #On gère le type avec une suite de if et on leur passe des paramètres
    return empty([], AbstractObserver)
end


struct TrajectoryMemoryDump <: AbstractDump
    n_save_iters::Int64
    time::Array
    xt::Array
    save_index::Int64
end

function TrajectoryMemoryDump(n_save_iters = 1; n_step = 1e4,kwargs...)
    n_save = floor(Int, n_step / n_save_iters)
    traj = Array{typeof(state.x[1])}(undef, size(state.x)[1])
    time = Array{Float64}(1:n_save)
    return TrajectoryMemoryDump(n_save_iters, time, traj, 1)
end

function run_obs(obs::TrajectoryMemoryDump, t::Float64, state::AbstractState; kwargs...)
    obs.xt[obs.save_index, :] .= deepcopy(state.x)
    obs.time[obs.save_index] = t
    obs.save_index += 1
end

struct TrajectoryFileDump <: AbstractDump
    n_save_iters::Int
    file_pattern::String
    file::IOStream
    function TrajectoryFileDump(n_save_iters,file_pattern; kwargs...)
        return TrajectoryFileDump(n_save_iters,file_pattern, stdout)
    end
end

function start_obs_traj(obs::TrajectoryFileDump)
    file = open(filename, "w")
end

function run_obs(obs::TrajectoryFileDump, t::Float64, state::AbstractState; kwargs...)
    writedlm(obs.file, vcat([t], state.x)')
end


function end_obs_traj(obs::TrajectoryFileDump)
    close(obs.file)
end


"""
Function that check if the stride is valid, if yes, perform the observation
"""
# function observation(obs::AbstractObserver,n::Int64,t::Float64,state::AbstractState)
#     if (mod(n, obs.n_save_iters) == 0)
#         run_obs(obs,t,state)
#     end
# end

# Defined to exist when there is no one
function run_obs(obs::AbstractObserver, t::Float64, state::AbstractState; kwargs...)

end

# Defined to exist when there is no one
function end_obs_traj(obs::AbstractObserver) # At the end of one traj

end

function close_obs(obs::AbstractObserver) # At the end of all traj

end
