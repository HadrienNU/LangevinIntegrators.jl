struct GJ{FP<:AbstractForce,TF<:AbstractFloat,TM} <: VelocityVerletIntegrator
    force::FP
    β::TF
    γ::TF
    M::TM
    Δt::TF
    sqrtM::TM
    a::TF
    b::TF
    σ::TF
    dim::Int64
    bc::Union{AbstractSpace,Nothing}
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


# Function spécifique pour GJF
function GJF(force::FP, β::TF, γ::TF, M::TM, Δt::TF,  dim::Int64=1, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
    a = (1 - 0.5 * γ * Δt) / (1 + 0.5 * γ * Δt)
    b = 1 / (1 + 0.5 * γ * Δt)
    σ = sqrt(2 * γ * Δt / β)
    sqrtM = sqrt.(M)
    return GJ(force, β, γ, M, Δt, sqrtM, a, b, σ, dim, bc)
end

#TODO: Implementer les autres formes de GJ
function GJ(force::FP, β::TF, γ::TF, M::TM, Δt::TF, type="I", dim::Int64=1, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
    #Faire un switch sur les valeur de type pour avoir les coeffs des autres GJ
    a = (1 - 0.5 * γ * Δt) / (1 + 0.5 * γ * Δt)
    b = 1 / (1 + 0.5 * γ * Δt)
    σ = sqrt(2 * γ * Δt / β)
    sqrtM = sqrt.(M)
    return GJ(force, β, γ, M, Δt, sqrtM, a, b, σ, dim, bc)
end

# A passer en forme compact
function UpdateState!(state::VelocityVerletState, integrator::GJ; kwargs...)

    state.ξ = randn(integrator.dim)

    state.x =
        state.x .+ integrator.b * integrator.Δt .* state.v .+ 0.5 * integrator.b * integrator.Δt^2 / integrator.M * state.f .+
        0.5 * integrator.b * integrator.Δt / integrator.sqrtM * integrator.σ * state.ξ
    apply_space!(integrator.bc,state.x,state.v)
    nostop = forceUpdate!(integrator.force, state.f_new, state.x; kwargs...)

    state.v =
        integrator.a * state.v .+ 0.5 * integrator.Δt / integrator.M * (integrator.a * state.f .+ state.f_new) .+
        integrator.b * integrator.sqrtM / integrator.M * integrator.σ * state.ξ

    state.f = state.f_new

    return nostop
end
