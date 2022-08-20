struct ABOBA_Hidden{TGV, TF<:AbstractFloat, TM} <: InertialIntegrator
    force::TGV
    Δt::TF
    S::AA
    friction_vv::Matrix{TF}
    friction_vh::Matrix{TF}
    friction_hv::Matrix{TF}
    friction_hh::Matrix{TF}
    dim::Int64
    dim_tot::Int64
end

"""
    ABOBA_Hidden(force, β, γ, M, Δt)

Set up the ABOBA_Hidden integrator for underdamped Langevin with hidden variables.

### Fields

* force   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
function ABOBA_Hidden(force::TGV, A::Array{TF},C::Array{TF}, Δt::TF,dim::Int)where{FP<:AbstractForce, TF<:AbstractFloat}
    dim_tot=size(A)
    friction = exp(-1 *  Δt * A)
    return ABOBA_Hidden(force, Δt, cholesky(friction*C+C*A').L,friction[1:dim,1:dim],friction[1:dim,1+dim:dim_tot],friction[1+dim:dim_tot,1:dim],friction[1+dim:dim_tot,1+dim:dim_tot],dim,dim_tot)
end

mutable struct HiddenABOBAState{TF<:AbstractFloat} <:AbstractMemoryHiddenState
    x::Vector{TF}
	v::Vector{TF}
    h::Vector{TF}
    q_mid::Vector{TF}
    p_mid::Vector{TF}
    p̂_mid::Vector{TF}
    f_mid::Vector{TF}
end

#TODO initialize velocity

function InitState!(x₀,v₀, integrator::ABOBA_Hidden)
    return HiddenABOBAState(x₀, v₀, similar(x₀), similar(x₀), similar(x₀), similar(x₀))
end

function InitState(x₀,v₀, integrator::ABOBA_Hidden)
    return HiddenABOBAState(deepcopy(x₀),deepcopy(v₀), similar(x₀), similar(x₀),similar(x₀), similar(x₀))
end

function InitState!(s::AbstractMemoryHiddenState, integrator::ABOBA_Hidden)
    return HiddenABOBAState(s.x,s.v, similar(s.x), similar(s.x),similar(s.x), similar(s.x))
end

function InitState(s::AbstractMemoryHiddenState, integrator::ABOBA_Hidden)
    return HiddenABOBAState(deepcopy(s.x),deepcopy(s.v), similar(s.x), similar(s.x),similar(s.x), similar(s.x))
end

function UpdateState!(state::HiddenABOBAState, integrator::ABOBA_Hidden)

    @. state.q_mid = state.x + 0.5 * integrator.Δt * state.v
    forceUpdate!(integrator.force,state.f_mid,state.q_mid)
    @. state.p_mid = state.v + 0.5 * integrator.Δt * state.f_mid

    gauss = integrator.S * randn(integrator.dim_tot) # For latter consider, putting gauss in state to reserve the memory

    @. state.p̂_mid = integrator.friction_vv* state.p_mid +integrator.friction_vh*state.h  + gauss[1:integrator.dim] # A remplacer par la ligne suivante
    @. state.h =  integrator.friction_hv*state.p_mid + integrator.friction_hh*state.h + gauss[1+integrator.dim:integrator.dim_tot]

    @. state.v = state.p̂_mid + 0.5 * integrator.Δt * state.f_mid
    @. state.x = state.q_mid + 0.5 * integrator.Δt * state.v




    state
end
