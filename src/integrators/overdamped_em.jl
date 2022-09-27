struct EM{FP<:AbstractForce,TF<:AbstractFloat} <: OverdampedIntegrator
    force::FP
    β::TF
    Δt::TF
    σ::TF
    bc::Union{AbstractSpace,Nothing}
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
function EM(force::FP, β::TF, Δt::TF, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat}
    σ = sqrt(2 * Δt / β)
    return EM(force, β, Δt, σ, bc)
end


mutable struct EMState{TF<:AbstractFloat} <: AbstractOverdampedState
    x::Vector{TF}
    f::Vector{TF}
    dim::Int64
    function EMState(x₀::Vector{TF}, f::Vector{TF}) where {TF<:AbstractFloat}
        return new{TF}(x₀, f, length(x₀))
    end
end


function InitState!(x₀, integrator::EM)
    f = forceUpdate(integrator.force, x₀)
    return EMState(x₀, f)
end

function InitState(x₀, integrator::EM)
    f = forceUpdate(integrator.force, x₀)
    return EMState(deepcopy(x₀), f)
end


function UpdateState!(state::EMState, integrator::EM; kwargs...)

    state.x = state.x .+ integrator.Δt .* state.f .+ integrator.σ .* randn(state.dim)
    apply_space!(integrator.bc,state.x)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)

    return nostop
end
