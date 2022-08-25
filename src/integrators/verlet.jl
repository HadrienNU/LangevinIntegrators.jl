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
end

mutable struct VerletState{TF<:AbstractFloat} <: AbstractInertialState
    x::Vector{TF}
    v::Vector{TF}
    f::Vector{TF}
    v_mid::Vector{TF}
    dim::Int64
end

function InitState!(x₀, v₀, integrator::Verlet)
    f = forceUpdate(integrator.force, x₀)
    return VerletState(x₀, v₀, f, similar(v₀), length(x₀))
end

function InitState(x₀, v₀, integrator::Verlet)
    f = forceUpdate(integrator.force, x₀)
    return VerletState(deepcopy(x₀), deepcopy(v₀), f, similar(v₀), length(x₀))
end

function UpdateState!(state::VerletState, integrator::Verlet; kwargs...)
    state.v_mid = state.v .+ 0.5 * integrator.Δt / integrator.M * state.f
    @. state.x = state.x + integrator.Δt * state.v_mid
    #apply_bc!(integrator.bc,state.x,state.v)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)
    state.v = state.v_mid .+ 0.5 * integrator.Δt / integrator.M * state.f
    return nostop
end
