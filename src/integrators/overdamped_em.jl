struct EM{FP<:AbstractForce, TF<:AbstractFloat} <: OverdampedIntegrator
    force::FP
    β::TF
    Δt::TF
    σ::TF
end
#Il faut faire aussi un EM_ND qui gère un sigma matriciel

"""
    EM(force, β, Δt)
Set up the EM integrator for overdamped Langevin.
### Fields
* force   - In place gradient of the potential
* β     - Inverse temperature
* Δt    - Time step
"""
function EM(force::FP, β::TF, Δt::TF) where{FP<:AbstractForce, TF<:AbstractFloat}
    σ = sqrt(2 * Δt /β)
    return EM(force, β, Δt, σ)
end


mutable struct EMState{TF<:AbstractFloat} <:AbstractOverdampedState
    x::Vector{TF}
    f::Vector{TF}
    dim::Int64
end


function InitState!(x₀, integrator::EM)
    f=forceUpdate(integrator.force, x₀)
    return EMState(x₀, copy(f),length(x₀))
end

function InitState(x₀, integrator::EM)
    f=forceUpdate(integrator.force, x₀)
    return EMState(deepcopy(x₀), copy(f),length(x₀))
end


function UpdateState!(state::EMState, integrator::EM)

    state.x = state.x + integrator.Δt * state.f + integrator.σ * randn(state.dim)
    # @timeit_debug timer "UpdateState: forceUpdate!" begin
        forceUpdate!(integrator.force,state.f, state.x)
    # end

    state
end
