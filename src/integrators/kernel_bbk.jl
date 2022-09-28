struct BBK_Kernel{FP<:AbstractForce,TF<:AbstractFloat,TM} <: KernelIntegrator
    force::FP
    β::TF
    kernel::Union{Vector{TF},Vector{Matrix{TF}}}
    σ_corr::Union{Vector{TF},Vector{Matrix{TF}}}
    M::TM
    Δt::TF
    dim::Int64
    bc::Union{AbstractSpace,Nothing}
end

"""
    BBK_Kernel(force, β, γ, M, Δt)

Set up the BBK_Kernel integrator for inertial Langevin.

Adapted from Iterative Reconstruction of Memory Kernels Gerhard Jung,*,†,‡ Martin Hanke,*,§ and Friederike Schmid*,†

### Fields

* force   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
function BBK_Kernel(force::FP, β::TF, kernel::Vector{Matrix{TF}}, M::TM, Δt::TF, dim::Int64=1, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
    noise_fdt=sqrt(Δt / β) / M * real.(ifft(sqrt.(fft(kernel)))) # note quand Kernel est une matrix il faut faire le cholesky

    return BBK_Kernel(force, β, kernel, noise_fdt, M, Δt, dim, bc)
end

mutable struct BBKKernelState{TF<:AbstractFloat} <: AbstractMemoryKernelState
    x::Vector{TF}
    v::Vector{TF}
    v_mid::Vector{TF}
    f::Vector{TF}
    diss_f::Vector{TF} # Dissipative force
    v_t::Queue{Vector{TF}} # Trajectory of v to compute the kernel
    noise_n::Queue{Vector{TF}} # Utiliser un Circular buffer à la place?

    function BBKKernelState(x₀::Vector{TF}, v₀::Vector{TF}, f::Vector{TF}, v_t, noise::Queue{Vector{TF}}) where {TF<:AbstractFloat}
        return new{TF}(x₀, v₀, similar(v₀), f, similar(f) , v_t , noise, length(x₀)) # TODO: import v_t as a vector and fill the queue by zeros and given v_t
    end
end

function InitState!(x₀, v₀, integrator::BBK_Kernel)
    f = forceUpdate(integrator.force, x₀)
    v_t=Queue{typeof(v₀)}()
    for n=1:length(integrator.kernel) # Check if this is the right size
        push!(v_t, zeros(state.dim))
    end
    noise=init_randn_correlated(length(integrator.σ_corr))
    return BBKKernelState(x₀, v₀, f, v_t, noise)
end

function InitState(x₀, v₀, integrator::BBK_Kernel)
    f = forceUpdate(integrator.force, x₀)
    v_t=Queue{typeof(v₀)}()
    for n=1:length(integrator.kernel) # Check if this is the right size
        push!(v_t, zeros(state.dim))
    end
    noise=init_randn_correlated(length(integrator.σ_corr))
    return BBKKernelState(deepcopy(x₀), deepcopy(v₀), f, v_t, noise)
end

function UpdateState!(state::BBKKernelState, integrator::BBK_Kernel; kwargs...)
    state.v_mid = state.v .+ 0.5 * integrator.Δt / integrator.M * state.f .-integrator.Δt .*( state.diss_f.+ 0.5*integrator.kernel[0]*state.v)
    @. state.x = state.x + integrator.Δt * state.v_mid
    apply_space!(integrator.bc,state.x,state.v)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)

    state.diss_f = corr_mem(state.v_t, integrator.kernel) .+ integrator.σ *  randn_correlated(state, integrator) # Pour ne calculer l'intégrale qu'une fois, on la stocke puisqu'elle resservira au prochain pas de temps
    state.v = (state.v_mid .+ 0.5 * integrator.Δt / integrator.M * state.f .+ integrator.Δt*state.diss_f) / (1 + 0.5 * integrator.Δt * integrator.kernel[0])

    return nostop
end
