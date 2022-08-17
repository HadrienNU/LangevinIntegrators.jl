struct BBK{TGV, TF<:AbstractFloat, TM} <: InertialIntegrator
    force::TGV
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
function BBK(force::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}
    σ = sqrt(γ * Δt / β)
    return BBK(force, β, γ, M, Δt, σ)
end

mutable struct BBKState{TF<:AbstractFloat} <:AbstractInertialState
    x::Vector{TF}
	v::Vector{TF}
    p_mid::Vector{TF}
    f::Vector{TF}
end

function InitState!(x₀, integrator::BBK)
    f=forceUpdate!(integrator.force, x₀[1])
    return BBKState(x₀, similar(x₀[1]), copy(f))
end

function InitState(x₀, integrator::BBK)

    f=forceUpdate!(integrator.force, x₀[1])
    return BBKState(deepcopy(x₀), similar(x₀[1]), copy(f))
end

function UpdateState!(state::BBKState, integrator::BBK)

    @. state.p_mid = state.v + 0.5 * integrator.Δt * state.f - 0.5 * integrator.Δt * integrator.γ/integrator.M * state.v + integrator.σ * randn()
    @. state.x = state.x + integrator.Δt/integrator.M * state.p_mid
    forceUpdate!(integrator.force,state.f, state.x)
    @. state.v = (state.p_mid + 0.5 * integrator.Δt * state.f + integrator.σ * randn())/(1 + 0.5 * integrator.Δt * integrator.γ/integrator.M)

    state
end
