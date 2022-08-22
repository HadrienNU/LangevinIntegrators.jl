struct ABOBA{FP<:AbstractForce, TF<:AbstractFloat, TM} <: InertialIntegrator
    force::FP
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
function ABOBA(force::FP, β::TF, γ::TF, M::TM, Δt::TF) where {FP<:AbstractForce, TF<:AbstractFloat,TM}

    c₀ = exp(-Δt * γ) / M
    c₁ = sqrt((1 - exp(-2*γ*Δt))/β)
    sqrtM = sqrt.(M) / M
    return ABOBA(force, β, γ, M, Δt, c₀, c₁, sqrtM)
end

mutable struct ABOBAState{TF<:AbstractFloat} <:AbstractInertialState
    x::Vector{TF}
	v::Vector{TF}
    x_mid::Vector{TF}
    v_mid::Vector{TF}
    p̂_mid::Vector{TF}
    f_mid::Vector{TF}
    dim::Int64
end

#TODO initialize velocity

function InitState!(x₀,v₀, integrator::ABOBA)
    return ABOBAState(x₀,v₀, similar(x₀), similar(v₀), similar(v₀), similar(x₀),length(x₀))
end

function InitState(x₀,v₀, integrator::ABOBA)
    return ABOBAState(deepcopy(x₀),deepcopy(v₀), similar(x₀), similar(v₀),similar(v₀), similar(x₀),length(x₀))
end

function UpdateState!(state::ABOBAState, integrator::ABOBA)

    @. state.x_mid = state.x + 0.5 * integrator.Δt * state.v
    nostop = forceUpdate!(integrator.force,state.f_mid,state.x_mid)
    state.v_mid = state.v .+ 0.5 * integrator.Δt/integrator.M * state.f_mid
    state.p̂_mid = integrator.c₀ .* state.v_mid + integrator.c₁ .* integrator.sqrtM * randn(state.dim)
    state.v = state.p̂_mid .+ 0.5 * integrator.Δt/integrator.M * state.f_mid
    @. state.x = state.x_mid + 0.5 * integrator.Δt * state.v

    return nostop
end
