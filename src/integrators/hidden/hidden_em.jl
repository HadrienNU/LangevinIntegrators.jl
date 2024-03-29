struct EM_Hidden{FP<:AbstractForce,TF<:AbstractFloat,AA<:AbstractArray} <: HiddenIntegrator
    force::FP
    Δt::TF
    σ::AA
    friction_vv::Array{TF}
    friction_vh::Array{TF}
    friction_hv::Array{TF}
    friction_hh::Array{TF}
    dim::Int64
    dim_tot::Int64
    bc::Union{AbstractSpace,Nothing}
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
function EM_Hidden(
    force::FP,
    A::Array{TF},
    C::Array{TF},
    Δt::TF,
    dim::Int,
    bc::Union{AbstractSpace,Nothing} = nothing,
) where {FP<:AbstractForce,TF<:AbstractFloat}
    dim_tot = size(A)[1]
    friction = A * Δt
    C_sym = 0.5 .* (C .+ C') #C should be symmetric
    return EM_Hidden(
        force,
        Δt,
        cholesky(friction * C_sym + C_sym * friction').L,
        friction[1:dim, 1:dim],
        friction[1:dim, (1+dim):dim_tot],
        friction[(1+dim):dim_tot, 1:dim],
        friction[(1+dim):dim_tot, (1+dim):dim_tot],
        dim,
        dim_tot,
        bc,
    )
end


mutable struct HiddenEMState{TF<:AbstractFloat} <: AbstractMemoryHiddenState
    x::Vector{TF}
    v::Vector{TF}
    h::Vector{TF}
    f::Vector{TF}
    friction_h::Vector{TF}
    ξ::Vector{TF}
    function HiddenEMState(
        x₀::Vector{TF},
        v₀::Vector{TF},
        h₀::Vector{TF},
        f::Vector{TF},
    ) where {TF<:AbstractFloat}
        return new{TF}(x₀, v₀, h₀, f, similar(h₀), similar(f))
    end
end

function InitState!(x₀, v₀, h₀, integrator::EM_Hidden)
    f = forceUpdate(integrator.force, x₀)
    return HiddenEMState(x₀, v₀, h₀, f)
end


function UpdateState!(state::HiddenEMState, integrator::EM_Hidden; kwargs...)

    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)
    state.x .+= integrator.Δt * state.v
    apply_space!(integrator.bc, state.x, state.v)
    state.ξ = integrator.σ * randn(integrator.dim_tot) # For latter consider, putting gauss in state to reserve the memory
    state.friction_h .=
        -integrator.friction_hv * state.v .- integrator.friction_hh * state.h
    @views state.v .+=
        -integrator.friction_vv * state.v .- integrator.friction_vh * state.h .+
        integrator.Δt * state.f .+ state.ξ[1:integrator.dim]
    @views state.h .+= state.friction_h .+ state.ξ[1+integrator.dim:integrator.dim_tot] #gauss and friction should be taking Dt into account

    return nostop
end
