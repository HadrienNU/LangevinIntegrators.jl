
"""
    VelocityVerlet(force, M, Δt)

Set up the velocity Verlet integrator.

### Fields

* force   - In place gradient of the potential
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
struct VelocityVerlet{
    FP<:AbstractForce,
    TF<:AbstractFloat,
    TFM<:Union{TF,AbstractMatrix{TF}},
} <: VelocityVerletIntegrator
    force::FP
    M::TFM
    Δt::TF
    σ::TFM
    dim::Int64
    bc::Union{AbstractSpace,Nothing}
    function VelocityVerlet(
        force::FP,
        M::TFM,
        Δt::TF,
        dim::Int64 = 1,
        bc::Union{AbstractSpace,Nothing} = nothing,
    ) where {FP<:AbstractForce,TF<:AbstractFloat,TFM<:Union{TF,AbstractMatrix{TF}}}
        new{FP,TF,TFM}(force, M, Δt, zero(M), dim, bc)
    end
end

function Verlet(
    force::FP,
    M::TFM,
    Δt::TF,
    dim::Int64 = 1,
    bc::Union{AbstractSpace,Nothing} = nothing,
) where {FP<:AbstractForce,TF<:AbstractFloat,TFM<:Union{TF,AbstractMatrix{TF}}}  # Equivalence of Verlet and VelocityVerlet
    return VelocityVerlet(force, M, Δt, dim, bc)
end


mutable struct VelocityVerletState{TF<:AbstractFloat} <: AbstractInertialState
    x::Vector{TF}
    v::Vector{TF}
    v_mid::Vector{TF}
    f::Vector{TF}
    ξ::Vector{TF}
    ξ₂::Vector{TF}
    function VelocityVerletState(
        x₀::Vector{TF},
        v₀::Vector{TF},
        f::Vector{TF},
        σ::TFM,
    ) where {TF<:AbstractFloat,TFM<:Union{TF,AbstractMatrix{TF}}}
        return new{TF}(x₀, v₀, similar(v₀), f, σ * randn(length(f)), σ * randn(length(f)))
    end
end


function InitState!(x₀, v₀, integrator::VVI) where {VVI<:VelocityVerletIntegrator}
    if integrator.dim != length(x₀)
        throw(ArgumentError("Mismatch of dimension in state initialization"))
    end
    f = forceUpdate(integrator.force, x₀)
    return VelocityVerletState(x₀, v₀, f, integrator.σ)
end


function UpdateState!(state::VelocityVerletState, integrator::VelocityVerlet; kwargs...)
    state.v_mid .= state.v .+ 0.5 * integrator.Δt / integrator.M * state.f
    @. state.x += integrator.Δt * state.v_mid
    apply_space!(integrator.bc, state.x, state.v_mid)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)
    state.v .= state.v_mid .+ 0.5 * integrator.Δt / integrator.M * state.f
    return nostop
end
