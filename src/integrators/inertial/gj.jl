struct GJ{FP<:AbstractForce,TF<:AbstractFloat,TM} <: VelocityVerletIntegrator
    force::FP
    β::TF
    γ::TF
    M::TM
    Δt::TF
    c2::TF
    sc1::TF
    sc3::TF
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


# Function spécifique pour GJF
function GJF(force::FP, β::TF, γ::TF, M::TM, Δt::TF,  dim::Int64=1, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
    return GJ(force, β, γ, M, Δt, "I", dim, bc)
end

#TODO: Implementer les autres formes de GJ
function GJ(force::FP, β::TF, γ::TF, M::TM, Δt::TF, type="I", dim::Int64=1, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
    #Faire un switch sur les valeur de type pour avoir les coeffs des autres GJ
    a = γ * Δt/ M
    if type == "I"
        c2 = (1 - 0.5* a) / (1 + 0.5 * a)
    elseif type == "II"
        c2 = exp(-a)
    elseif type == "III"
        c2 = 1-a
    elseif type == "IV"
        c2 = (sqrt(1+4*a)-1)/(2*a)
    elseif type == "V"
        c2 = 1/(1+a)
    elseif type == "VI"
        c2 = 1/(1+0.5*a)^2
    else # Raise an error
        println("Unknown GJ type")
    end
    sc1 = sqrt((1+c2)/2)
    sc3  = sqrt((1-c2)/a)
    σ = sqrt(2 * γ * Δt / β) / M
    return GJ(force, β, γ, M, Δt, c2, sc1, sc3, σ, dim, bc)
end

# A passer en forme compact
function UpdateState!(state::VelocityVerletState, integrator::GJ; kwargs...)

    state.ξ = integrator.σ *randn(integrator.dim)

    state.v_mid = integrator.sc1 * state.v + 0.5*integrator.sc3 *integrator.Δt*state.f/integrator.M + 0.5*integrator.sc3*state.ξ

    state.x = state.x .+ integrator.sc3*integrator.Δt * state.v_mid

    apply_space!(integrator.bc,state.x,state.v)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)

    state.v = integrator.c2/integrator.sc1 * state.v_mid +  0.5*(integrator.sc3/integrator.sc1 )* integrator.Δt *state.f/integrator.M  + 0.5*(integrator.sc3/integrator.sc1)*state.ξ



    return nostop
end
