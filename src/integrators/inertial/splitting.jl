struct BAOAB{FP<:AbstractForce,TF<:AbstractFloat,TFM<:Union{TF,AbstractMatrix{TF}}, NP <: AbstractNoise} <:
       VelocityVerletIntegrator
    Δt::TF
    force::FP
    M::Union{TF,TFM}
    β::TF
    γ::TFM
    noise::NP
    c₂::TFM
    σ::TFM
    dim::Int64
    bc::Union{AbstractSpace,Nothing}
end

"""
    BAOAB(force, β, γ, M, Δt)

Set up the BAOAB integrator for inertial Langevin.

### Fields

* force   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
function BAOAB(
    force::FP,
    β::TF,
    γ::Union{TF,TM},
    M::Union{TF,TM},
    Δt::TF,
    dim::Int64 = 1,
    noise= nothing:: Union{Nothing,AbstractNoise},
    bc::Union{AbstractSpace,Nothing} = nothing,
) where {FP<:AbstractForce,TF<:AbstractFloat,TM<:AbstractMatrix{TF}}
    c₂ = exp(-Δt * γ) / M
    σ = sqrt((1 - exp(-2 * γ * Δt)) / β) / sqrt(M)
    if isnothing(noise)
        noise=GaussianNoise(dim)
    end
    return BAOAB(Δt,force, M, β, γ, noise, c₂, σ, dim, bc)
end


struct OBABO{FP<:AbstractForce,TF<:AbstractFloat,TFM<:Union{TF,AbstractMatrix{TF}}, NP <: AbstractNoise} <:
       VelocityVerletIntegrator
    Δt::TF
    force::FP
    M::Union{TF,TFM}
    β::TF
    γ::TFM
    noise::NP
    noise2::NP
    cc₂::TFM
    c₂::TFM
    σ::TFM
    dim::Int64
    bc::Union{AbstractSpace,Nothing}
end


"""
    OBABO(force, β, γ, M, Δt)

Set up the OBABO integrator for inertial Langevin.

### Fields

* force   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
function OBABO(
    force::FP,
    β::TF,
    γ::Union{TF,TM},
    M::Union{TF,TM},
    Δt::TF,
    dim::Int64 = 1,
    noise= nothing:: Union{Nothing,AbstractNoise},
    bc::Union{AbstractSpace,Nothing} = nothing,
) where {FP<:AbstractForce,TF<:AbstractFloat,TM<:AbstractMatrix{TF}}
    cc₂ = exp(-0.5 * Δt * γ) / M
    c₂ = exp(-Δt * γ) / M
    σ = sqrt((1 - exp(-γ * Δt)) / β) / sqrt(M)
    if isnothing(noise)
        noise=GaussianNoise(dim)
        # Et sion faire des test pour savoir si tout est cohérent
    end  # Attention il doit y avoir 2 bruits pour OBABO
    return OBABO(Δt,force, M, β, γ, noise, deepcopy(noise), cc₂, c₂, σ, dim, bc)
end



function UpdateState!(state::VelocityVerletState, integrator::BAOAB; kwargs...)

    state.v_mid .= state.v .+ 0.5 * integrator.Δt / integrator.M * state.f
    generate_noise!(state.ξ, integrator.noise, state.x, state.v_mid)
    state.ξ .= integrator.σ * state.ξ
    # state.ξ .= integrator.σ * randn(integrator.dim)
    state.x .+=
        0.5 * integrator.Δt * state.v_mid .+
        0.5 * integrator.Δt * (integrator.c₂ .* state.v_mid .+ state.ξ)
    apply_space!(integrator.bc, state.x, state.v)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)
    state.v .=
        integrator.c₂ * state.v_mid .+ 0.5 * integrator.Δt / integrator.M * state.f .+
        state.ξ

    return nostop
end


# TODO: Ajouter 2 bruits
function UpdateState!(state::VelocityVerletState, integrator::OBABO; kwargs...)
    generate_noise!(state.ξ₂, integrator.noise2, state.x, state.v_mid)
    state.ξ₂ .= integrator.σ * state.ξ₂
    # state.ξ₂ .= integrator.σ * randn(integrator.dim)
    state.v_mid .=
        integrator.cc₂ * state.v .+ 0.5 * integrator.Δt / integrator.M * state.f .+ state.ξ₂
    @. state.x += integrator.Δt * state.v_mid
    apply_space!(integrator.bc, state.x, state.v)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)
    generate_noise!(state.ξ, integrator.noise, state.x, state.v_mid)
    state.ξ .= integrator.σ * state.ξ
    # state.ξ .= integrator.σ * randn(integrator.dim)
    state.v .=
        integrator.cc₂ * (state.v_mid .+ 0.5 * integrator.Δt / integrator.M * state.f) .+
        state.ξ

    return nostop
end
