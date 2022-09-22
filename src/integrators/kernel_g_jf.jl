struct Kernel_GJF{FP<:AbstractForce,TF<:AbstractFloat,TM} <: KernelIntegrator
    force::FP
    β::TF
    kernel::Union{Vector{TF},Vector{Matrix{TF}}}
    σ_corr::Union{Vector{TF},Vector{Matrix{TF}}}
    M::TM
    Δt::TF
    sqrtM::TM
    a::Union{TF,Matrix{TF}}
    b::Union{TF,Matrix{TF}}
    bc::Union{AbstractSpace,Nothing}
end

"""
    Kernel_GJF(force, β, γ, M, Δt)

Set up the G-JF integrator for inertial Langevin.

Adapted from Iterative Reconstruction of Memory Kernels Gerhard Jung,*,†,‡ Martin Hanke,*,§ and Friederike Schmid*,†

### Fields

* force   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
function Kernel_GJF(force::FP, β::TF, γ::TF, M::TM, Δt::TF, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
    a = (1 - 0.5 * γ * Δt) / (1 + 0.5 * γ * Δt)
    b = 1 / (1 + 0.5 * γ * Δt)
    σ = sqrt(2 * γ * Δt / β)
    sqrtM = sqrt.(M)
    return Kernel_GJF(force, β, γ, M, Δt, sqrtM, a, b, σ, bc)
end

mutable struct GJFKernelState{TF<:AbstractFloat} <: AbstractMemoryKernelState
    x::Vector{TF}
    v::Vector{TF}
    f::Vector{TF}
    f_new::Vector{TF}
    x_t::Vector{Vector{TF}} # Trajectory of x to compute the kernel, both array are given by the size of the kernel
    W_t::Vector{Vector{TF}} # Utiliser un Circular buffer à la place?
    dim::Int64
end

function InitState!(x₀, v₀, integrator::Kernel_GJF)
    f = forceUpdate(integrator.force, x₀)
    return GJFKernelState(x₀, v₀, f, copy(f), similar(f), length(x₀))
end

function InitState(x₀, v₀, integrator::Kernel_GJF)
    f = forceUpdate(integrator.force, x₀)
    return GJFKernelState(deepcopy(x₀), deepcopy(v₀), f, copy(f), similar(f), length(x₀))
end

function UpdateState!(state::GJFKernelState, integrator::Kernel_GJF; kwargs...)

    state.ξ = randn(state.dim) # To replace with correlated noise

    mem_int = corr_mem(state.x_t[1:] - state.x_t[:-1] , integrator.kernel) # Check limit to give, les  limites sont pas bonnes pour x_t

    state.x =
        state.x .+ integrator.b * integrator.Δt .* state.v .+ 0.5 * integrator.b * integrator.Δt^2 / integrator.M * state.f
        .- 0.5 * integrator.b * integrator.Δt * mem_int .+ 0.5 * integrator.b * integrator.Δt / integrator.sqrtM * integrator.σ * state.ξ

    apply_space!(integrator.bc,state.x,state.v)
    nostop = forceUpdate!(integrator.force, state.f_new, state.x; kwargs...)

    state.v =
        integrator.a * state.v .+ 0.5 * integrator.Δt / integrator.M * (integrator.a * state.f .+ state.f_new)
        .- integrator.b * mem_int .+ integrator.b * integrator.sqrtM / integrator.M * integrator.σ * state.ξ

    state.f = state.f_new

    return nostop
end
