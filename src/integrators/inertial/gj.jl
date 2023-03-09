struct GJ{FP<:AbstractForce,TF<:AbstractFloat,TM} <: VelocityVerletIntegrator
    force::FP
    β::TF
    γ::TF
    M::TM
    Δt::TF
    c₂::TF
    sc₁::TF
    d₁::TF
    σ::TF
    dim::Int64
    bc::Union{AbstractSpace,Nothing}
end

"""
    GJF(force, β, γ, M, Δt)

Set up the G-JF integrator for inertial Langevin.

### Fields

* force   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
function GJF(force::FP, β::TF, γ::TF, M::TM, Δt::TF,  dim::Int64=1, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
    return GJ(force, β, γ, M, Δt, "I", dim, bc)
end

"""
    GJ(force, β, γ, M, Δt, type)

Set up the various GJ integrator for inertial Langevin.
Use type to select the approriate integrator

### Fields

* force   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
* type - Choice of the integrator should be one of "I","II","III","IV","V","VI"
"""
function GJ(force::FP, β::TF, γ::TF, M::TM, Δt::TF, type="I", dim::Int64=1, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
    #Faire un switch sur les valeur de type pour avoir les coeffs des autres GJ
    a = γ * Δt/ M
    if type == "I"
        c₂ = (1 - 0.5* a) / (1 + 0.5 * a)
    elseif type == "II"
        c₂ = exp(-a)
    elseif type == "III"
        c₂ = 1-a
    elseif type == "IV"
        c₂ = (sqrt(1+4*a)-1)/(2*a)
    elseif type == "V"
        c₂ = 1/(1+a)
    elseif type == "VI"
        c₂ = 1/(1+0.5*a)^2
    else # Raise an error
        println("Unknown GJ type")
    end
    sc₁ = sqrt((1+c₂)/2)
    d₁  = sqrt((1-c₂)/a)
    σ = sqrt(2 * γ * Δt / β) / sqrt(M)
    return GJ(force, β, γ, M, Δt, c₂, sc₁, d₁, σ, dim, bc)
end


function UpdateState!(state::VelocityVerletState, integrator::GJ; kwargs...)

    state.ξ = integrator.σ *randn(integrator.dim)

    state.v_mid = integrator.sc₁ * state.v + 0.5*integrator.d₁ *integrator.Δt*state.f/integrator.M + 0.5*integrator.d₁*state.ξ

    state.x = state.x .+ integrator.d₁*integrator.Δt * state.v_mid

    apply_space!(integrator.bc,state.x,state.v)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)

    state.v = (integrator.c₂ * state.v_mid +  0.5*integrator.d₁* integrator.Δt *state.f/integrator.M  + 0.5*integrator.d₁*state.ξ)/integrator.sc₁



    return nostop
end
