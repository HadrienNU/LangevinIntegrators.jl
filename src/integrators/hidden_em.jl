struct EM_Hidden{FP<:AbstractForce, TF<:AbstractFloat, AA <: AbstractArray} <: HiddenIntegrator
    force::FP
    Δt::TF
    S::AA
    friction_vv::Array{TF}
    friction_vh::Array{TF}
    friction_hv::Array{TF}
    friction_hh::Array{TF}
    dim::Int64
    dim_tot::Int64
end

"""
    EM_Hidden(force, β, Δt)
Set up the EM_Hidden integrator for underdamped Langevin with hidden variables.
### Fields
* force   - In place gradient of the potential
* A     - Friction matrix
* C     - diffusion matrix
* Δt    - Time step
"""
#Evidemment à changer
function EM_Hidden(force::FP, A::Array{TF},C::Array{TF}, Δt::TF,dim::Int) where{FP<:AbstractForce, TF<:AbstractFloat}
    dim_tot=size(A)
    friction = A * Δt
    return EM_Hidden(force, Δt, cholesky(friction*C+C*A').L,friction[1:dim,1:dim],friction[1:dim,1+dim:dim_tot],friction[1+dim:dim_tot,1:dim],friction[1+dim:dim_tot,1+dim:dim_tot],dim,dim_tot)
end


mutable struct HiddenEMState{TF<:AbstractFloat} <:AbstractMemoryHiddenState
    x::Vector{TF}
    v::Vector{TF}
    h::Vector{TF}
    f::Vector{TF}

    # friction_h::Vector{TF}
end

function InitState!(x₀, v₀, h₀, integrator::EM_Hidden)
    f=forceUpdate(integrator.force, x₀)
    return HiddenEMState(x₀, v₀, h₀, copy(f))
end

function InitState(x₀, v₀, h₀, integrator::EM_Hidden)
    f=forceUpdate(integrator.force, x₀)
    return HiddenEMState(deepcopy(x₀), deepcopy(v₀), deepcopy(h₀),copy(f))
end

function InitState!(s::AbstractMemoryHiddenState, integrator::EM_Hidden)
    f=forceUpdate(integrator.force, s.x)
    return HiddenEMState(s.x,s.v,s.h, copy(f))
end

function InitState(s::AbstractMemoryHiddenState, integrator::EM_Hidden)
    f=forceUpdate(integrator.force, s.x)
    return HiddenEMState(deepcopy(s.x), deepcopy(s.v),deepcopy(s.h),copy(f))
end

function UpdateState!(state::HiddenEMState, integrator::EM_Hidden)

    state.x = state.x + integrator.Δt * state.v / integrator.M
    forceUpdate!(integrator.force, state.f, state.x)

    gauss = integrator.S * randn(integrator.dim_tot) # For latter consider, putting gauss in state to reserve the memory
    friction_h = - integrator.friction_hv*state.v - integrator.friction_hh*state.h
    state.v = state.v - integrator.friction_vv* state.v -integrator.friction_vh*state.h + integrator.Δt * state.f + gauss[1:integrator.dim]
    state.h = state.h  + friction_h + gauss[1+integrator.dim:integrator.dim_tot] #gauss and friction should be taking Dt into account

    state
end
