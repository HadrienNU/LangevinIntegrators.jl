struct BAOAB{FP<:AbstractForce,TF<:AbstractFloat,TM} <: VelocityVerletIntegrator
    force::FP
    β::TF
    γ::TF
    M::TM
    Δt::TF
    c₀::TF
    c₁::TF
    sqrtM::TM
    dim::Int64
    bc::Union{AbstractSpace,Nothing}
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
function BAOAB(force::FP, β::TF, γ::TF, M::TM, Δt::TF, dim::Int64=1, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
    c₀ = exp(-Δt * γ) / M
    c₁ = sqrt((1 - exp(-2 * γ * Δt)) / β)
    sqrtM = sqrt.(M) / M
    return BAOAB(force, β, γ, M, Δt, c₀, c₁, sqrtM, dim, bc)
end


struct OBABO{FP<:AbstractForce,TF<:AbstractFloat,TM} <: VelocityVerletIntegrator
    force::FP
    β::TF
    γ::TF
    M::TM
    Δt::TF
    c₀::TF
    c₁::TF
    sqrtM::TM
    dim::Int64
    bc::Union{AbstractSpace,Nothing}
end


"""
    OBABO(force, β, γ, M, Δt)

Set up the OBABO integrator for inertial Langevin.

### Fields

* force   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
#TODO: A écrire
function OBABO(force::FP, β::TF, γ::TF, M::TM, Δt::TF, dim::Int64=1, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
    c₀ = exp(-Δt * γ) / M
    c₁ = sqrt((1 - exp(-2 * γ * Δt)) / β)
    sqrtM = sqrt.(M) / M
    return OBABO(force, β, γ, M, Δt, c₀, c₁, sqrtM, dim, bc)
end



function UpdateState!(state::VelocityVerletState, integrator::BAOAB; kwargs...)

    state.v_mid = state.v .+ 0.5 * integrator.Δt / integrator.M * state.f
    state.ξ = randn(integrator.dim)
    @. state.x = state.x .+ 0.5 * integrator.Δt * state.v_mid + 0.5 * integrator.Δt * (integrator.c₀ .* state.v_mid .+ integrator.c₁ .* integrator.sqrtM * state.ξ)
    apply_space!(integrator.bc,state.x,state.v)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)
    state.v = integrator.c₀ * state.v_mid  + 0.5 * integrator.Δt / integrator.M * state.f .+ integrator.c₁ * integrator.sqrtM * state.ξ

    return nostop
end


# A écrire là c'est pour BAOAB
function UpdateState!(state::VelocityVerletState, integrator::OBABO; kwargs...)

    state.v_mid = state.v .+ 0.5 * integrator.Δt / integrator.M * state.f
    state.ξ = randn(integrator.dim)
    @. state.x = state.x .+ 0.5 * integrator.Δt * state.v_mid + 0.5 * integrator.Δt * (integrator.c₀ .* state.v_mid .+ integrator.c₁ .* integrator.sqrtM * state.ξ)
    apply_space!(integrator.bc,state.x,state.v)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)
    state.v = integrator.c₀ * state.v_mid  + 0.5 * integrator.Δt / integrator.M * state.f .+ integrator.c₁ * integrator.sqrtM * state.ξ

    return nostop
end
