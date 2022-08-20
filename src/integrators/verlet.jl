"""
    Verlet(force, M, Δt)

Set up the Verlet integrator.

### Fields

* force   - In place gradient of the potential
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
struct Verlet{TGV, TF<:AbstractFloat, TM} <: InertialIntegrator
    force::TGV
    M::TM
    Δt::TF
end

mutable struct VerletState{TF<:AbstractFloat} <:AbstractInertialState
    x::Vector{TF}
	v::Vector{TF}
    f::Vector{TF}
    v_mid::Vector{TF}
end

function InitState!(x₀, integrator::Verlet)
    f= forceUpdate(integrator.force,x₀)
    return VerletState(x₀, copy(f) , similar(x₀[1]));
end

function InitState(x₀, integrator::Verlet)
	f= forceUpdate(integrator.force,x₀)
    return VerletState(deepcopy(x₀), copy(f) , similar(x₀[1]));
end

function UpdateState!(state::VerletState, integrator::Verlet)
    @. state.v_mid = state.v + 0.5 * integrator.Δt * state.f/integrator.M;
    @. state.x = state.x + integrator.Δt * state.v_mid;
	forceUpdate!(integrator.force,state.f,state.x)
    @. state.v = state.v_mid + 0.5 * integrator.Δt * state.f/integrator.M;
    state
end
