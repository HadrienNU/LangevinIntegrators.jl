struct BAOAB{TGV, TF<:AbstractFloat, TM} <: InertialIntegrator
    force::TGV
    β::TF
    γ::TF
    M::TM
    Δt::TF
    c₀::TF
    c₁::TF
    sqrtM::TM
end

"""
    BAOAB(force, β, γ, M, Δt)

Set up the BAOAB integrator for inertial Langevin.

### Fields

* force   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
function BAOAB(force::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}

    c₀ = exp(-Δt * γ)
    c₁ = sqrt((1 - exp(-2*γ*Δt))/β)
    sqrtM = sqrt.(M)
    return BAOAB(force, β, γ, M, Δt, c₀, c₁, sqrtM)
end

mutable struct BAOABState{TF<:AbstractFloat} <:AbstractInertialState
    x::Vector{TF}
	v::Vector{TF}
    q_mid::Vector{TF}
    p_mid::Vector{TF}
    p̂_mid::Vector{TF}
    f::Vector{TF}
end

function InitState!(x₀, integrator::BAOAB)
    f=forceUpdate!(integrator.force, x₀[1])
    return BAOABState(x₀, similar(x₀[1]), similar(x₀[1]),similar(x₀[1]), copy(f))
end

function InitState(x₀, integrator::BAOAB)

    f=forceUpdate!(integrator.force, x₀[1])
    return BAOABState(deepcopy(x₀), similar(x₀[1]), similar(x₀[1]),similar(x₀[1]), copy(f))
end

function UpdateState!(state::BAOABState, integrator::BAOAB)

    @. state.p_mid = state.v - 0.5 * integrator.Δt * state.f
    @. state.q_mid = state.x + 0.5 * integrator.Δt * state.p_mid/integrator.M
    @. state.p̂_mid = integrator.c₀ * state.p_mid + integrator.c₁ * integrator.sqrtM * randn()
    @. state.x = state.q_mid + 0.5 * integrator.Δt * state.p̂_mid/integrator.M
    forceUpdate!(integrator.force,state.f,state.x)
    @. state.v = state.p̂_mid - 0.5 * integrator.Δt * state.f

    state
end
