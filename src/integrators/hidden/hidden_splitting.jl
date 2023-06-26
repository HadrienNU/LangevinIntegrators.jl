struct BAOAB_Hidden{FP<:AbstractForce,TF<:AbstractFloat,TM} <:
       HiddenVelocityVerletIntegrator
    force::FP
    C::Matrix{TF}
    γ::Matrix{TF}
    M::TM
    Δt::TF
    c₂::Array{TF}
    c₂_vh::Array{TF}
    c₂_hv::Array{TF}
    c₂_hh::Array{TF}
    σ::Array{TF}
    dim::Int64
    dim_tot::Int64
    bc::Union{AbstractSpace,Nothing}
end

"""
    BAOAB_Hidden(force, β, γ, M, Δt)

Set up the BAOAB_Hidden integrator for underdamped Langevin with hidden variables.

### Fields

* force   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
function BAOAB_Hidden(
    force::FP,
    A::Array{TF},
    C::Array{TF},
    M::TM,
    Δt::TF,
    dim::Int,
    bc::Union{AbstractSpace,Nothing} = nothing,
) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
    dim_tot = size(A)[1]
    friction = exp(-1 * Δt * A)
    C_sym = 0.5 .* (C .+ C')
    return BAOAB_Hidden(
        force,
        C,
        A,
        M,
        Δt,
        friction[1:dim, 1:dim],
        friction[1:dim, (1+dim):dim_tot],
        friction[(1+dim):dim_tot, 1:dim],
        friction[(1+dim):dim_tot, (1+dim):dim_tot],
        cholesky(friction * C_sym + C_sym * friction').L,
        dim,
        dim_tot,
        bc,
    )
end

mutable struct HiddenVelocityVerletState{TF<:AbstractFloat} <: AbstractMemoryHiddenState
    x::Vector{TF}
    v::Vector{TF}
    h::Vector{TF}
    h_mid::Vector{TF}
    v_mid::Vector{TF}
    v̂_mid::Vector{TF}
    f::Vector{TF}
    ξ::Vector{TF}
    ξ₂::Vector{TF}
    function HiddenVelocityVerletState(
        x₀::Vector{TF},
        v₀::Vector{TF},
        h₀::Vector{TF},
        f::Vector{TF},
    ) where {TF<:AbstractFloat}
        return new{TF}(x₀, v₀, h₀, similar(x₀), similar(v₀), similar(v₀), f)
    end
end

function InitState!(x₀, v₀, h₀, integrator::HVI) where {HVI<:HiddenVelocityVerletIntegrator}
    if integrator.dim != length(x₀)
        throw(ArgumentError("Mismatch of dimension in state initialization"))
    end
    f = forceUpdate(integrator.force, x₀)
    return HiddenVelocityVerletState(x₀, v₀, h₀, f, integrator.σ)
end

function UpdateState!(state::VelocityVerletState, integrator::BAOAB_Hidden; kwargs...)

    state.v_mid .= state.v .+ 0.5 * integrator.Δt / integrator.M * state.f
    state.ξ .= integrator.σ * randn(integrator.dim_tot)

    state.v̂_mid =
        integrator.c₂ .* state.v_mid .+ integrator.c₂_vh * state.h .+
        state.ξ[1:integrator.dim]
    state.h =
        integrator.c₂_hv * state.v_mid .+ integrator.c₂_hh * state.h .+
        state.ξ[1+integrator.dim:integrator.dim_tot]
    state.x .+= 0.5 * integrator.Δt * state.v_mid .+ 0.5 * integrator.Δt * state.v̂_mid
    apply_space!(integrator.bc, state.x, state.v)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)

    state.v .= state.v̂_mid .+ 0.5 * integrator.Δt / integrator.M * state.f


    return nostop
end
