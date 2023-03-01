struct BBK{FP<:AbstractForce,TF<:AbstractFloat,TM} <: VelocityVerletIntegrator
    force::FP
    β::TF
    γ::TF
    M::TM
    Δt::TF
    c₀::TF
    c₁::TF
    c₂::TF
    σ::TF
    dim::Int64
    bc::Union{AbstractSpace,Nothing}
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
function BBK(force::FP, β::TF, γ::TF, M::TM, Δt::TF, dim::Int64=1, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
    c₀ = (1 - 0.5 * Δt * γ / M)
    c₁ = 1.0 / (1 + 0.5 * Δt * γ / M)
    σ = sqrt(2 * γ * Δt / β) / sqrt(M)
    c₂ = c₀*c₁
    return BBK(force, β, γ, M, Δt, c₀, c₁, c₂, σ, dim, bc)
end

struct VEC{FP<:AbstractForce,TF<:AbstractFloat,TM} <: VelocityVerletIntegrator
    force::FP
    β::TF
    γ::TF
    M::TM
    Δt::TF
    c₁::TF
    c₂::TF
    d₁::TF
    d₂::TF
    σ::TF
    dim::Int64
    bc::Union{AbstractSpace,Nothing}
end

"""
    VEC(force, β, γ, M, Δt)

Set up the Vanden-Eijnden Ciccotti integrator for inertial Langevin.
Taken from "Second-order integrators for Langevin equations with holonomic constraints" doi: 10.1016/j.cplett.2006.07.086

### Fields

* force   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
function VEC(force::FP, β::TF, γ::TF, M::TM, Δt::TF, dim::Int64=1, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
    c₂ = (1 - 0.5 * γ * Δt/ M +  0.125*(γ * Δt)^2/ M)
    c₁ =  0.5 * Δt *(1 - 0.25 * γ * Δt) / M
    d₁ = 0.5 *(1 - 0.25 * γ * Δt)
    d₂ = - 0.25 * γ * Δt/sqrt(3)
    σ = sqrt(2 * γ * Δt / β) / sqrt(M)
    return VEC(force, β, γ, M, Δt, c₁, c₂, d₁, d₂, σ, dim, bc)
end


function UpdateState!(state::VelocityVerletState, integrator::BBK; kwargs...)
    # state.ξ = integrator.σ * randn(integrator.dim)
    state.v_mid = integrator.c₀ * state.v .+ 0.5 * integrator.Δt / integrator.M * state.f  .+ 0.5*state.ξ
    @. state.x = state.x + integrator.Δt * state.v_mid
    apply_space!(integrator.bc,state.x,state.v)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)
    state.ξ = integrator.σ * randn(integrator.dim)
    state.v = integrator.c₁ *(state.v_mid .+ 0.5 * integrator.Δt / integrator.M * state.f + 0.5*state.ξ)
    return nostop
end




function UpdateState!(state::VelocityVerletState, integrator::VEC; kwargs...)
    state.ξ = integrator.σ * randn(integrator.dim)
    state.ξ₂ = integrator.σ * randn(integrator.dim)
    state.v_mid = integrator.c₂ * state.v .+ integrator.c₁* state.f  .+ integrator.d₁*state.ξ + integrator.d₂*state.ξ₂
    @. state.x = state.x + integrator.Δt * (state.v_mid + (0.5/sqrt(3))*state.ξ₂)
    apply_space!(integrator.bc,state.x,state.v)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)
    state.v = integrator.c₂ * state.v_mid .+ integrator.c₁* state.f  .+ integrator.d₁*state.ξ + integrator.d₂*state.ξ₂
    return nostop
end
