struct EM_Hidden{FP<:AbstractForce, TF<:AbstractFloat} <: HiddenIntegrator
    force::FP
    β::TF
    Δt::TF
    σ::TF
end
#Il faut faire aussi un EM_ND qui gère un sigma matriciel

"""
    EM(force, β, Δt)
Set up the EM integrator for overdamped Langevin.
### Fields
* force   - In place gradient of the potential
* β     - Inverse temperature
* Δt    - Time step
"""
function EM_Hidden(force::FP, β::TF, Δt::TF) where{FP<:AbstractForce, TF<:AbstractFloat}
    σ = sqrt(2 * Δt /β)
    return EM_Hidden(force, β, Δt, σ)
end

#Pour le EM_1D et EM_ND, vérifier que la dim dans la force est la même que celle passé

#Ce qu'on peut faire c'est une struct EMState1D et EMStateND pour différencier le cas 1D qui est le plus fréquent
#Et alors on a juste à écrire un force update différents et surcharger la function similar pour les floats

abstract type  AbstractEMHiddenState <: AbstractOverdampedState end

mutable struct EMState{TF<:AbstractFloat} <:AbstractEMHiddenState
    x::Vector{TF}
    f::Vector{TF}
end


function InitState!(x₀, integrator::EM)
    f=forceUpdate(integrator.force, x₀)
    return EMState(x₀, copy(f))
end

function InitState(x₀, integrator::EM)
    f=forceUpdate(integrator.force, x₀)
    return EMState(deepcopy(x₀), copy(f))
end

function InitState!(s::AbstractOverdampedState, integrator::EM)
    f=forceUpdate(integrator.force, s.x)
    return EMState(s.x, copy(f))
end

function InitState(s::AbstractOverdampedState, integrator::EM)
    f=forceUpdate(integrator.force, s.x)
    return EMState(deepcopy(s.x), copy(f))
end

function UpdateState!(state::EMState, integrator::EM)

    @. state.x = state.x + integrator.Δt * state.f + integrator.σ * randn()
    forceUpdate!(integrator.force,state.f, state.x)

    state
end
