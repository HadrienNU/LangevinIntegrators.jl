struct EM_Kernel{FP<:AbstractForce,TF<:AbstractFloat,AA<:AbstractArray} <: KernelIntegrator
    force::FP
    Δt::TF
    kernel::Union{Vector{TF},Vector{Matrix{TF}}}
    σ_corr::Union{Vector{TF},Vector{Matrix{TF}}}
    dim::Int64
    bc::Union{AbstractSpace,Nothing}
end

"""
    EM_Kernel(force, β, Δt)
Set up the EM_Kernel integrator for underdamped Langevin with hidden variables.
### Fields
* force   - In place gradient of the potential
* A     - Friction matrix
* C     - diffusion matrix
* Δt    - Time step
"""
#Evidemment à changer
function EM_Kernel(force::FP,  β::TF, kernel::Union{Vector{TF},Vector{Matrix{TF}}}, Δt::TF, dim::Int64=1, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat}
    noise_fdt=sqrt(Δt / β) * real.(ifft(sqrt.(fft(kernel)))) # note quand Kernel est une matrix il faut faire le cholesky
    #σ_corr, doit inclure le Δt
    return EM_Kernel(force, Δt, kernel, noise_fdt, dim,  bc)
end


mutable struct KernelEMState{TF<:AbstractFloat} <: AbstractMemoryKernelState
    x::Vector{TF}
    v::Vector{TF}
    f::Vector{TF}
    v_t::Queue{Vector{TF}} # Trajectory of v to compute the kernel
    noise_n::Queue{Vector{TF}}
    dim::Int64
    # friction_h::Vector{TF}
    function KernelEMState(x₀::Vector{TF}, v₀::Vector{TF}, f::Vector{TF}, v_t, noise::Queue{Vector{TF}}) where {TF<:AbstractFloat}
        return new{TF}(x₀, v₀, f, v_t , noise, length(x₀))
    end
end

function InitState!(x₀, v₀, integrator::EM_Kernel)
    f = forceUpdate(integrator.force, x₀)
    v_t=Queue{typeof(v₀)}()
    for n=1:length(integrator.kernel) # Check if this is the right size
        push!(v_t, zeros(state.dim))
    end
    noise=init_randn_correlated(length(integrator.σ_corr))
    return KernelEMState(x₀, v₀, h₀, f, v_t, noise)
end

function InitState(x₀, v₀, integrator::EM_Kernel)
    f = forceUpdate(integrator.force, x₀)
    v_t=Queue{typeof(v₀)}()
    for n=1:length(integrator.kernel) # Check if this is the right size
        push!(v_t, zeros(state.dim))
    end
    noise=init_randn_correlated(length(integrator.σ_corr))
    return KernelEMState(deepcopy(x₀), deepcopy(v₀), f, v_t)
end



function UpdateState!(state::KernelEMState, integrator::EM_Kernel; kwargs...)
    mem_int= corr_mem(state.v_t,integrator.kernel)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)
    @. state.x = state.x + integrator.Δt * state.v
    apply_space!(integrator.bc,state.x,state.v)
    state.v = state.v .+integrator.Δt / integrator.M * state.f .+ integrator.Δt*mem_ker .+ randn_correlated(state,integrator)
    return nostop
end
