struct BAOAB{FP<:AbstractForce, TF<:AbstractFloat, TM} <: InertialIntegrator
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
    BAOAB(force, β, γ, M, Δt)

Set up the BAOAB integrator for inertial Langevin.

### Fields

* force   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
function BAOAB(force::FP, β::TF, γ::TF, M::TM, Δt::TF) where {FP<:AbstractForce, TF<:AbstractFloat,TM}

    c₀ = exp(-Δt * γ) / M
    c₁ = sqrt((1 - exp(-2*γ*Δt))/β)
    sqrtM = sqrt.(M) / M
    return BAOAB(force, β, γ, M, Δt, c₀, c₁, sqrtM)
end

mutable struct BAOABState{TF<:AbstractFloat} <:AbstractInertialState
    x::Vector{TF}
	v::Vector{TF}
    x_mid::Vector{TF}
    v_mid::Vector{TF}
    p̂_mid::Vector{TF}
    f::Vector{TF}
    dim::Int64
end

function InitState!(x₀,v₀, integrator::BAOAB)
    f=forceUpdate(integrator.force, x₀)
    return BAOABState(x₀,v₀, similar(x₀), similar(v₀),similar(v₀), f,length(x₀))
end

function InitState(x₀,v₀, integrator::BAOAB)
    f=forceUpdate(integrator.force, x₀)
    return BAOABState(deepcopy(x₀),deepcopy(v₀), similar(x₀), similar(x₀),similar(x₀), f,length(x₀))
end

function UpdateState!(state::BAOABState, integrator::BAOAB)

    state.v_mid = state.v .+ 0.5 * integrator.Δt /integrator.M * state.f
    @. state.x_mid = state.x .+ 0.5 * integrator.Δt * state.v_mid
    #apply_bc!(integrator.bc,state.x_mid,state.v)
    state.p̂_mid = integrator.c₀ .* state.v_mid .+ integrator.c₁ .* integrator.sqrtM * randn(state.dim)
    @. state.x = state.x_mid + 0.5 * integrator.Δt * state.p̂_mid
    #apply_bc!(integrator.bc,state.x,state.v)
    nostop = forceUpdate!(integrator.force,state.f,state.x)
    state.v = state.p̂_mid .+ 0.5 * integrator.Δt/integrator.M * state.f

    return nostop
end
