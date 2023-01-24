struct BBK{FP<:AbstractForce,TF<:AbstractFloat,TM} <: VelocityVerletIntegrator
    force::FP
    β::TF
    γ::TF
    M::TM
    Δt::TF
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
    σ = sqrt(γ * Δt / β) / M
    return BBK(force, β, γ, M, Δt, σ, dim, bc)
end


function UpdateState!(state::VelocityVerletState, integrator::BBK; kwargs...)
    # state.ξ = integrator.σ * randn(integrator.dim)
    state.v_mid = state.v .+ 0.5 * integrator.Δt / integrator.M * state.f .- 0.5 * integrator.Δt .* integrator.γ * state.v .+ state.ξ
    @. state.x = state.x + integrator.Δt * state.v_mid
    apply_space!(integrator.bc,state.x,state.v)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)
    state.ξ = integrator.σ * randn(integrator.dim)
    state.v = (state.v_mid .+ 0.5 * integrator.Δt / integrator.M * state.f + state.ξ) / (1 + 0.5 * integrator.Δt * integrator.γ)
    return nostop
end
