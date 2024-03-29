struct BAOAB{FP<:AbstractForce,TF<:AbstractFloat,TM} <: VelocityVerletIntegrator
    force::FP
    β::TF
    γ::TF
    M::TM
    Δt::TF
    c₂::TF
    σ::TF
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
function BAOAB(
    force::FP,
    β::TF,
    γ::TF,
    M::TM,
    Δt::TF,
    dim::Int64 = 1,
    bc::Union{AbstractSpace,Nothing} = nothing,
) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
    c₂ = exp(-Δt * γ) / M
    σ = sqrt((1 - exp(-2 * γ * Δt)) / β) / sqrt(M)
    return BAOAB(force, β, γ, M, Δt, c₂, σ, dim, bc)
end


struct OBABO{FP<:AbstractForce,TF<:AbstractFloat,TM} <: VelocityVerletIntegrator
    force::FP
    β::TF
    γ::TF
    M::TM
    Δt::TF
    cc₂::TF
    c₂::TF
    σ::TF
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
function OBABO(
    force::FP,
    β::TF,
    γ::TF,
    M::TM,
    Δt::TF,
    dim::Int64 = 1,
    bc::Union{AbstractSpace,Nothing} = nothing,
) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
    cc₂ = exp(-0.5 * Δt * γ) / M
    c₂ = exp(-Δt * γ) / M
    σ = sqrt((1 - exp(-γ * Δt)) / β) / sqrt(M)
    return OBABO(force, β, γ, M, Δt, cc₂, c₂, σ, dim, bc)
end



function UpdateState!(state::VelocityVerletState, integrator::BAOAB; kwargs...)

    state.v_mid .= state.v .+ 0.5 * integrator.Δt / integrator.M * state.f
    state.ξ .= integrator.σ * randn(integrator.dim)
    state.x .+=
        0.5 * integrator.Δt * state.v_mid .+
        0.5 * integrator.Δt * (integrator.c₂ .* state.v_mid .+ state.ξ)
    apply_space!(integrator.bc, state.x, state.v)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)
    state.v .=
        integrator.c₂ * state.v_mid .+ 0.5 * integrator.Δt / integrator.M * state.f .+
        state.ξ

    return nostop
end



function UpdateState!(state::VelocityVerletState, integrator::OBABO; kwargs...)

    state.ξ₂ .= integrator.σ * randn(integrator.dim)
    state.v_mid .=
        integrator.cc₂ * state.v .+ 0.5 * integrator.Δt / integrator.M * state.f .+ state.ξ₂
    @. state.x += integrator.Δt * state.v_mid
    apply_space!(integrator.bc, state.x, state.v)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)
    state.ξ .= integrator.σ * randn(integrator.dim)
    state.v .=
        integrator.cc₂ * (state.v_mid .+ 0.5 * integrator.Δt / integrator.M * state.f) .+
        state.ξ

    return nostop
end
