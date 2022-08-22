struct BBK{FP<:AbstractForce, TF<:AbstractFloat, TM} <: InertialIntegrator
    force::FP
    β::TF
    γ::TF
    M::TM
    Δt::TF
    σ::TF
end

"""
    BBK(force, β, γ, M, Δt)

Set up the BBK integrator for inertial Langevin.

### Fields

* force   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
function BBK(force::FP, β::TF, γ::TF, M::TM, Δt::TF) where {FP<:AbstractForce, TF<:AbstractFloat,TM}
    σ = sqrt(γ * Δt / β) / M
    return BBK(force, β, γ, M, Δt, σ)
end

mutable struct BBKState{TF<:AbstractFloat} <:AbstractInertialState
    x::Vector{TF}
	v::Vector{TF}
    v_mid::Vector{TF}
    f::Vector{TF}
    dim::Int64
end

function InitState!(x₀, integrator::BBK)
    f=forceUpdate!(integrator.force, x₀[1])
    return BBKState(x₀, similar(x₀[1]), copy(f),length(x₀))
end

function InitState(x₀, integrator::BBK)

    f=forceUpdate!(integrator.force, x₀[1])
    return BBKState(deepcopy(x₀), similar(x₀[1]), copy(f),length(x₀))
end

function UpdateState!(state::BBKState, integrator::BBK)

    state.v_mid = state.v + 0.5 * integrator.Δt / integrator.M * state.f - 0.5 * integrator.Δt * integrator.γ * state.v + integrator.σ * randn(state.dim)
    state.x = state.x + integrator.Δt * state.v_mid
    forceUpdate!(integrator.force,state.f, state.x)
    state.v = (state.v_mid + 0.5 * integrator.Δt / integrator.M * state.f + integrator.σ * randn(state.dim))/(1 + 0.5 * integrator.Δt * integrator.γ)

    state
end
