struct ABOBA{FP<:AbstractForce,TF<:AbstractFloat,TM} <: InertialIntegrator
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
    ABOBA(force, β, γ, M, Δt)

Set up the ABOBA integrator for inertial Langevin.

### Fields

* force   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
function ABOBA(force::FP, β::TF, γ::TF, M::TM, Δt::TF, dim::Int64=1, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat,TM}

    c₀ = exp(-Δt * γ) / M
    c₁ = sqrt((1 - exp(-2 * γ * Δt)) / β)
    sqrtM = sqrt.(M) / M
    return ABOBA(force, β, γ, M, Δt, c₀, c₁, sqrtM, dim, bc)
end

mutable struct ABOBAState{TF<:AbstractFloat} <: AbstractInertialState
    x::Vector{TF}
    v::Vector{TF}
    x_mid::Vector{TF}
    v_mid::Vector{TF}
    p̂_mid::Vector{TF}
    f_mid::Vector{TF}
    function ABOBAState(x₀::Vector{TF}, v₀::Vector{TF}, f::Vector{TF}) where {TF<:AbstractFloat}
        return new{TF}(x₀, v₀, similar(x₀), similar(v₀), similar(v₀), f)
    end
end

#TODO initialize velocity

function InitState!(x₀, v₀, integrator::ABOBA)
    if integrator.dim != length(x₀)
        throw(ArgumentError("Mismatch of dimension in state initialization"))
    end
    f = forceUpdate(integrator.force, x₀)
    return ABOBAState(x₀, v₀, f)
end

function UpdateState!(state::ABOBAState, integrator::ABOBA; kwargs...)

    @. state.x_mid = state.x + 0.5 * integrator.Δt * state.v
    apply_space!(integrator.bc,state.x_mid,state.v)
    nostop = forceUpdate!(integrator.force, state.f_mid, state.x_mid; kwargs...)
    # Au passage il faudra rajouter ici un terme de métrique quand il est présent
    state.v_mid = state.v .+ 0.5 * integrator.Δt / integrator.M * state.f_mid
    state.p̂_mid = integrator.c₀ .* state.v_mid + integrator.c₁ .* integrator.sqrtM * randn(integrator.dim)
    state.v = state.p̂_mid .+ 0.5 * integrator.Δt / integrator.M * state.f_mid
    @. state.x = state.x_mid + 0.5 * integrator.Δt * state.v
    apply_space!(integrator.bc,state.x,state.v)

    return nostop
end
