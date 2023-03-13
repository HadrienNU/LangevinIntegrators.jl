
"""
    PositionVerlet(force, M, Δt)

Set up the position Verlet integrator.

### Fields

* force   - In place gradient of the potential
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
struct PositionVerlet{FP<:AbstractForce,TF<:AbstractFloat, TFM <: Union{TF,AbstractMatrix{TF}}} <: PositionVerletIntegrator
    force::FP
    M::TFM
    Δt::TF
    σ::TFM
    dim::Int64
    bc::Union{AbstractSpace,Nothing}
    function PositionVerlet(force::FP, M::TFM, Δt::TF, dim::Int64=1, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat, TFM <: Union{TF,AbstractMatrix{TF}}}
        new{FP,TF,TFM}(force, M, Δt, zero(M), dim, bc)
    end
end


struct ABOBA{FP<:AbstractForce,TF<:AbstractFloat, TFM <: Union{TF,AbstractMatrix{TF}}} <: PositionVerletIntegrator
    force::FP
    β::TF
    γ::TFM
    M::Union{TF,TFM}
    Δt::TF
    c₂::TFM
    σ::TFM
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
function ABOBA(force::FP, β::TF, γ::Union{TF,TM}, M::Union{TF,TM}, Δt::TF, dim::Int64=1, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat, TM <: AbstractMatrix{TF}}
    c₂ = exp(-Δt * γ) / M
    σ = sqrt((1 - exp(-2 * γ * Δt)) / β) / sqrt(M)
    return ABOBA(force, β, γ, M, Δt, c₂, σ, dim, bc)
end

mutable struct PositionVerletState{TF<:AbstractFloat} <: AbstractInertialState
    x::Vector{TF}
    v::Vector{TF}
    x_mid::Vector{TF}
    f_mid::Vector{TF}
    ξ::Vector{TF}
    function PositionVerletState(x₀::Vector{TF}, v₀::Vector{TF}, f::Vector{TF}, σ::TFM) where {TF<:AbstractFloat, TFM <: Union{TF,AbstractMatrix{TF}}}
        return new{TF}(x₀, v₀, similar(x₀), f, σ*randn(length(f)))
    end
end


function InitState!(x₀, v₀, integrator::PVI ) where {PVI <:PositionVerletIntegrator}
    if integrator.dim != length(x₀)
        throw(ArgumentError("Mismatch of dimension in state initialization"))
    end
    f = forceUpdate(integrator.force, x₀)
    return PositionVerletState(x₀, v₀, f, integrator.σ)
end

function UpdateState!(state::PositionVerletState, integrator::PositionVerlet; kwargs...)

    @. state.x_mid = state.x + 0.5 * integrator.Δt * state.v
    apply_space!(integrator.bc,state.x_mid,state.v)
    nostop = forceUpdate!(integrator.force, state.f_mid, state.x_mid; kwargs...)
    # Au passage il faudra rajouter ici un terme de métrique quand il est présent
    state.v .+= integrator.Δt / integrator.M * state.f_mid
    @. state.x = state.x_mid + 0.5 * integrator.Δt * state.v
    apply_space!(integrator.bc,state.x,state.v)
    return nostop
end


function UpdateState!(state::PositionVerletState, integrator::ABOBA; kwargs...)

    @. state.x_mid = state.x + 0.5 * integrator.Δt * state.v
    apply_space!(integrator.bc,state.x_mid,state.v)
    nostop = forceUpdate!(integrator.force, state.f_mid, state.x_mid; kwargs...)
    # Au passage il faudra rajouter ici un terme de métrique quand il est présent
    state.ξ .= integrator.σ *randn(integrator.dim)
    state.v .=  integrator.c₂ * (state.v .+ 0.5 * integrator.Δt / integrator.M * state.f_mid) .+ 0.5 * integrator.Δt / integrator.M * state.f_mid .+  state.ξ
    @. state.x = state.x_mid + 0.5 * integrator.Δt * state.v
    apply_space!(integrator.bc,state.x,state.v)

    return nostop
end
