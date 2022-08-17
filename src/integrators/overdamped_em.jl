struct EM{FP<:AbstractForce, TF<:AbstractFloat} <: OverdampedIntegrator
    force::FP
    β::TF
    Δt::TF
    σ::TF
end

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
end


function InitState!(x₀, integrator::EM)
    f=forceUpdate(integrator.force, x₀)
    return EMState(x₀, copy(f))
end

function InitState(x₀, integrator::EM)
    f=forceUpdate(integrator.force, x₀)
    return EMState(deepcopy(x₀), copy(f))
end

function InitState!(s::AbstractOverdampedState, integrator::EM)
    f=forceUpdate(integrator.force, s.x)
    return EMState(s.x, copy(f))
end

function InitState(s::AbstractOverdampedState, integrator::EM)
    f=forceUpdate(integrator.force, s.x)
    return EMState(deepcopy(s.x), copy(f))
end

function UpdateState!(state::EMState, integrator::EM)

    @. state.x = state.x + integrator.Δt * state.f + integrator.σ * randn()
    forceUpdate!(integrator.force,state.f, state.x)

    state
end
