struct GJF{TGV, TF<:AbstractFloat, TM} <: InertialIntegrator
    force::TGV
    β::TF
    γ::TF
    M::TM
    Δt::TF
    sqrtM::TM
    a::TF
    b::TF
    σ::TF
end

"""
    GJF(force, β, γ, M, Δt)

Set up the G-JF integrator for inertial Langevin.

### Fields

* force   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
function GJF(force::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}
    a = (1 - 0.5 * γ * Δt)/(1 + 0.5 * γ * Δt)
    b = 1 / (1 + 0.5 * γ * Δt)
    σ = sqrt(2 * γ * Δt / β)
    sqrtM = sqrt.(M)
    return GJF(force, β, γ, M, Δt, sqrtM, a, b, σ)
end

mutable struct GJFState{TF<:AbstractFloat} <:AbstractInertialState
    x::Vector{TF}
	v::Vector{TF}
    f::Vector{TF}
    f_new::Vector{TF}
    ξ::Vector{TF}
end

function InitState!(x₀, integrator::GJF)
    f = similar(x₀)
    f=forceUpdate!(integrator.force, x₀)
    return GJFState(x₀, similar(x₀), copy(f), similar(x₀), similar(x₀))
end

function InitState(x₀, integrator::GJF)

    f = similar(x₀)
    f=forceUpdate!(integrator.force, x₀)
    return GJFState(deepcopy(x₀), similar(x₀), copy(f), similar(x₀), similar(x₀))
end

function UpdateState!(state::GJFState, integrator::GJF)

    @. state.ξ = randn()

    @. state.x = state.x + integrator.b * integrator.Δt * state.v + 0.5 * integrator.b * integrator.Δt^2 / integrator.M * state.f + 0.5 * integrator.b * integrator.Δt / integrator.sqrtM * integrator.σ * state.ξ

    forceUpdate!(integrator.force,state.f_new, state.x)

    @. state.v = integrator.a * state.v + 0.5 * integrator.Δt/ integrator.M  * (integrator.a * state.f + state.f_new) + integrator.b * integrator.sqrtM / integrator.M  * integrator.σ * state.ξ

    @. state.f = state.f_new

    state
end
