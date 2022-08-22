#=
Les functions pour le dump
=#

abstract type AbstractObserver end

abstract type AbstractDump <: AbstractObserver end
abstract type AbstractStatisticalObs <: AbstractObserver end

struct TrajectoryMemoryDump <: AbstractDump
    n_save_iters::Int64
    time::Array
    xt::Array
    save_index::Int64
end

function TrajectoryMemoryDump(params::LangevinParams)
    n_save_iters=params.n_save_iters
    n_save=floor(Int, params.n_iters / n_save_iters)
    traj = Array{typeof(state.x[1])}(undef,(,size(state.x)[1]))
    time= Array{Float64}(1:n_save)
    return TrajectoryMemoryDump(n_save_iters,time,traj,1)
end

function run_obs(obs::TrajectoryMemoryDump,t::Float64,state::AbstractState;kwargs...)
    obs.xt[obs.save_index,:] .= deepcopy(state.x)
    obs.time[obs.save_index] = t
    obs.save_index += 1
end

struct TrajectoryFileDump <: AbstractDump
    n_save_iters::Int
    file::IOStream
end

function TrajectoryFileDump(params::LangevinParams)
    n_save_iters=params.n_save_iters
    file=open(filename,"w")
    return TrajectoryFileDump(n_save_iters,file)
end

function run_obs(obs::TrajectoryFileDump,t::Float64,state::AbstractState;kwargs...)
    writedlm(obs.file, vcat([t],state.x)')
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
function run_obs(obs::AbstractObserver,t::Float64,state::AbstractState;kwargs...)

end

# Defined to exist when there is no one
function end_obs_traj(obs::AbstractObserver)

end
