"""
    Verlet(force, M, Δt)

Set up the Verlet integrator.

### Fields

* force   - In place gradient of the potential
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
struct Verlet{FP<:AbstractForce,TF<:AbstractFloat,TM} <: InertialIntegrator
    force::FP
    M::TM
    Δt::TF
    dim::Int64
    bc::Union{AbstractSpace,Nothing}
    function Verlet(force::FP, M::TM, Δt::TF, dim::Int64=1, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
        new{FP,TF,TM}(force, M, Δt, dim, bc)
    end
end




mutable struct VerletState{TF<:AbstractFloat} <: AbstractInertialState
    x::Vector{TF}
    v::Vector{TF}
    f::Vector{TF}
    v_mid::Vector{TF}
    function VerletState(x₀::Vector{TF}, v₀::Vector{TF}, f::Vector{TF}) where {TF<:AbstractFloat}
        return new{TF}(x₀, v₀, f, similar(v₀))
    end
end

function InitState!(x₀, v₀, integrator::Verlet)
    f = forceUpdate(integrator.force, x₀)
    #Add check on dimension
    if integrator.dim != length(x₀)
        throw(ArgumentError("Mismatch of dimension in state initialization"))
    end
    return VerletState(x₀, v₀, f)
end

function UpdateState!(state::VerletState, integrator::Verlet; kwargs...)
    state.v_mid = state.v .+ 0.5 * integrator.Δt / integrator.M * state.f
    @. state.x = state.x + integrator.Δt * state.v_mid
    apply_space!(integrator.bc,state.x, state.v_mid)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)
    state.v = state.v_mid .+ 0.5 * integrator.Δt / integrator.M * state.f
    return nostop
end
