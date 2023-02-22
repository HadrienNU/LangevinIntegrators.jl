struct ABOBA_Hidden{FP<:AbstractForce,TF<:AbstractFloat,AA<:AbstractArray} <: HiddenIntegrator
    force::FP
    Δt::TF
    S::AA
    friction_vv::Matrix{TF}
    friction_vh::Matrix{TF}
    friction_hv::Matrix{TF}
    friction_hh::Matrix{TF}
    dim::Int64
    dim_tot::Int64
    bc::Union{AbstractSpace,Nothing}
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
function ABOBA_Hidden(force::FP, A::Array{TF}, C::Array{TF}, Δt::TF, dim::Int, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat}
    dim_tot = size(A)[1]
    friction = exp(-1 * Δt * A)
    C_sym = 0.5 .* (C .+ C')
    return ABOBA_Hidden(
        force,
        Δt,
        cholesky(friction * C_sym + C_sym * friction').L,
        friction[1:dim, 1:dim],
        friction[1:dim, (1+dim):dim_tot],
        friction[(1+dim):dim_tot, 1:dim],
        friction[(1+dim):dim_tot, (1+dim):dim_tot],
        dim,
        dim_tot,
        bc
    )
end

mutable struct HiddenABOBAState{TF<:AbstractFloat} <: AbstractMemoryHiddenState
    x::Vector{TF}
    v::Vector{TF}
    h::Vector{TF}
    x_mid::Vector{TF}
    p_mid::Vector{TF}
    p̂_mid::Vector{TF}
    f_mid::Vector{TF}
    function HiddenABOBAState(x₀::Vector{TF}, v₀::Vector{TF}, h₀::Vector{TF},f::Vector{TF}) where {TF<:AbstractFloat}
        return new{TF}(x₀, v₀, h₀, similar(x₀), similar(v₀), similar(v₀), f)
    end
end

#TODO initialize velocity

function InitState!(x₀, v₀, h₀, integrator::ABOBA_Hidden)
    f = forceUpdate(integrator.force, x₀)
    return HiddenABOBAState(x₀, v₀, h₀, f)
end

function UpdateState!(state::HiddenABOBAState, integrator::ABOBA_Hidden; kwargs...)

    @. state.x_mid = state.x + 0.5 * integrator.Δt * state.v
    apply_space!(integrator.bc,state.x_mid,state.v)
    nostop = forceUpdate!(integrator.force, state.f_mid, state.x_mid; kwargs...)
    @. state.p_mid = state.v + 0.5 * integrator.Δt * state.f_mid

    gauss = integrator.S * randn(integrator.dim_tot) # For latter consider, putting gauss in state to reserve the memory

    state.p̂_mid = integrator.friction_vv * state.p_mid .+ integrator.friction_vh * state.h .+ gauss[1:integrator.dim] # A remplacer par la ligne suivante
    state.h = integrator.friction_hv * state.p_mid .+ integrator.friction_hh * state.h .+ gauss[1+integrator.dim:integrator.dim_tot]

    @. state.v = state.p̂_mid + 0.5 * integrator.Δt * state.f_mid
    @. state.x = state.x_mid + 0.5 * integrator.Δt * state.v
    apply_space!(integrator.bc,state.x,state.v)

    return nostop
end
