struct EM_Kernel{FP<:AbstractForce,TF<:AbstractFloat,AA<:AbstractArray} <: KernelIntegrator
    force::FP
    Δt::TF
    kernel::Union{Vector{TF},Vector{Matrix{TF}}}
    σ_corr::Union{Vector{TF},Vector{Matrix{TF}}}
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
function EM_Kernel(force::FP, kernel::Union{Vector{TF},Vector{Matrix{TF}}}, Δt::TF, dim::Int, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat}

    #σ_corr, doit inclure le Δt
    return EM_Kernel(
        force,
        Δt,
        kernel,
        noise_fdt,
        bc
    )
end


mutable struct KernelEMState{TF<:AbstractFloat} <: AbstractMemoryKernelState
    x::Vector{TF}
    v::Vector{TF}
    h::Vector{TF}
    f::Vector{TF}
    v_t::Vector{Vector{TF}} # Trajectory of v to compute the kernel # Utiliser un Circular buffer à la place?
    noise_t::Vector{Vector{TF}} # Utiliser un Circular buffer à la place?
    dim::Int64
    # friction_h::Vector{TF}
end

function InitState!(x₀, v₀, h₀, integrator::EM_Kernel)
    f = forceUpdate(integrator.force, x₀)
    return KernelEMState(x₀, v₀, h₀, f, length(x₀))
end

function InitState(x₀, v₀, h₀, integrator::EM_Kernel)
    f = forceUpdate(integrator.force, x₀)
    return KernelEMState(deepcopy(x₀), deepcopy(v₀), deepcopy(h₀), f, length(x₀))
end



function UpdateState!(state::KernelEMState, integrator::EM_Kernel; kwargs...)
    mem_int= corr_mem(state.v_t,integrator.kernel)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)
    @. state.x = state.x + integrator.Δt * state.v
    apply_space!(integrator.bc,state.x,state.v)
    state.v = state.v .+integrator.Δt / integrator.M * state.f .+ integrator.Δt*mem_ker .+ randn_correlated(state,integrator)
    return nostop
end
