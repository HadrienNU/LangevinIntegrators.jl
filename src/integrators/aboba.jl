struct ABOBA{TGV, TF<:AbstractFloat, TM} <: InertialIntegrator
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
    ABOBA(force, β, γ, M, Δt)

Set up the ABOBA integrator for inertial Langevin.

### Fields

* force   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
function ABOBA(force::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}

    c₀ = exp(-Δt * γ)
    c₁ = sqrt((1 - exp(-2*γ*Δt))/β)
    sqrtM = sqrt.(M)
    return ABOBA(force, β, γ, M, Δt, c₀, c₁, sqrtM)
end

mutable struct ABOBAState{TF<:AbstractFloat} <:AbstractInertialState
    x::Vector{TF}
	v::Vector{TF}
    q_mid::Vector{TF}
    p_mid::Vector{TF}
    p̂_mid::Vector{TF}
    f_mid::Vector{TF}
end

#TODO initialize velocity

function InitState!(x₀,v₀, integrator::ABOBA)
    return ABOBAState(x₀,v₀, similar(x₀), similar(x₀), similar(x₀), similar(x₀))
end

function InitState(x₀,v₀, integrator::ABOBA)
    return ABOBAState(deepcopy(x₀),deepcopy(v₀), similar(x₀), similar(x₀),similar(x₀), similar(x₀))
end

function InitState!(s::AbstractInertialState, integrator::EM)
    return ABOBAState(s.x,s.v, similar(s.x), similar(s.x),similar(s.x), similar(s.x))
end

function InitState(s::AbstractInertialState, integrator::EM)
    return ABOBAState(deepcopy(s.x),deepcopy(s.v), similar(s.x), similar(s.x),similar(s.x), similar(s.x))
end

function UpdateState!(state::ABOBAState, integrator::ABOBA)

    @. state.q_mid = state.x + 0.5 * integrator.Δt * state.v/integrator.M
    forceUpdate!(integrator.force,state.f_mid,state.q_mid)
    @. state.p_mid = state.v + 0.5 * integrator.Δt * state.f_mid
    @. state.p̂_mid = integrator.c₀ * state.p_mid + integrator.c₁ * integrator.sqrtM * randn()
    @. state.v = state.p̂_mid + 0.5 * integrator.Δt * state.f_mid
    @. state.x = state.q_mid + 0.5 * integrator.Δt * state.v/integrator.M

    state
end
