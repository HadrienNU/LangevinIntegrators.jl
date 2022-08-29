struct GJF{FP<:AbstractForce,TF<:AbstractFloat,TM} <: InertialIntegrator
    force::FP
    β::TF
    γ::TF
    M::TM
    Δt::TF
    sqrtM::TM
    a::TF
    b::TF
    σ::TF
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
function GJF(force::FP, β::TF, γ::TF, M::TM, Δt::TF, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
    a = (1 - 0.5 * γ * Δt) / (1 + 0.5 * γ * Δt)
    b = 1 / (1 + 0.5 * γ * Δt)
    σ = sqrt(2 * γ * Δt / β)
    sqrtM = sqrt.(M)
    return GJF(force, β, γ, M, Δt, sqrtM, a, b, σ, bc)
end

mutable struct GJFState{TF<:AbstractFloat} <: AbstractInertialState
    x::Vector{TF}
    v::Vector{TF}
    f::Vector{TF}
    f_new::Vector{TF}
    ξ::Vector{TF}
    dim::Int64
end

function InitState!(x₀, v₀, integrator::GJF)
    f = forceUpdate(integrator.force, x₀)
    return GJFState(x₀, v₀, f, copy(f), similar(f), length(x₀))
end

function InitState(x₀, v₀, integrator::GJF)
    f = forceUpdate(integrator.force, x₀)
    return GJFState(deepcopy(x₀), deepcopy(v₀), f, copy(f), similar(f), length(x₀))
end

function UpdateState!(state::GJFState, integrator::GJF; kwargs...)

    state.ξ = randn(state.dim)

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
