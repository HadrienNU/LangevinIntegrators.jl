struct EM{FP<:AbstractForce, TF<:AbstractFloat} <: OverdampedIntegrator
    force::FP
    β::TF
    Δt::TF
    σ::TF
    # ou alors de BC global comme les contraintes, on va coder les 2 et avec les types paramétriques ça va marcher
    # On a un type BC_indep qui est égal à Array{BC} et un type Contraintes pour le premier c'est element-wise
    # bc::Array{BC}  # On a une BC par dimension, es-ce que ça ne serait pas une propriétésde l'intégrateur
    #Si les BC sont dans l'intégrateur, on a pas besoin d'actualiser la force
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
function EM(force::FP, β::TF, Δt::TF) where{FP<:AbstractForce, TF<:AbstractFloat}
    σ = sqrt(2 * Δt /β)
    return EM(force, β, Δt, σ)
end


mutable struct EMState{TF<:AbstractFloat} <:AbstractOverdampedState
    x::Vector{TF}
    f::Vector{TF}
    dim::Int64
end


function InitState!(x₀, integrator::EM)
    f=forceUpdate(integrator.force, x₀)
    return EMState(x₀, copy(f),length(x₀))
end

function InitState(x₀, integrator::EM)
    f=forceUpdate(integrator.force, x₀)
    return EMState(deepcopy(x₀), copy(f),length(x₀))
end


function UpdateState!(state::EMState, integrator::EM)

    state.x = state.x .+ integrator.Δt .* state.f .+ integrator.σ .* randn(state.dim)
    #apply_bc!(integrator.bc,state.x)
    # @timeit_debug timer "UpdateState: forceUpdate!" begin
    nostop = forceUpdate!(integrator.force,state.f, state.x)
    # end

    return nostop
end
