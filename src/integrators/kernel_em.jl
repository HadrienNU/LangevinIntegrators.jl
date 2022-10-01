struct EM_Kernel{FP<:AbstractForce,TF<:AbstractFloat,AA<:AbstractArray} <: KernelIntegrator
    force::FP
    Δt::TF
    kernel::Array{TF}
    σ_corr::Array{TF}
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
    if  size(kernel,2) != dim || size(kernel,2) != size(kernel,3)
        throw(ArgumentError("Mismatch of dimension bewteen kernel and space dimension"))
    end
    ker_mat=reshape(kernel,dim,dim, :)
    noise_fdt=sqrt(Δt / β) * real.(ifft(sqrt.(fft(ker_mat,3)))) # note quand Kernel est une matrix il faut faire le cholesky
    #σ_corr, doit inclure le Δt
    return EM_Kernel(force, Δt, ker_mat, noise_fdt, dim,  bc)
end


mutable struct KernelEMState{TF<:AbstractFloat} <: AbstractMemoryKernelState
    x::Vector{TF}
    v::Vector{TF}
    f::Vector{TF}
    v_t::Vector{Vector{TF}} # Trajectory of v to compute the kernel
    noise_n::Vector{Vector{TF}}
    # friction_h::Vector{TF}
    function KernelEMState(x₀::Vector{TF}, v₀::Vector{TF}, f::Vector{TF}, v_t, noise::Vector{Vector{TF}}) where {TF<:AbstractFloat}
        return new{TF}(x₀, v₀, f, v_t , noise)
    end
end

function InitState!(x₀, v₀, integrator::EM_Kernel)
    f = forceUpdate(integrator.force, x₀)
    v_t=Vector{typeof(v₀)}()
    push!(v_t, v₀)
    noise=init_randn_correlated(integrator.σ_corr)
    return KernelEMState(x₀, v₀, h₀, f, v_t, noise)
end

function InitState(x₀, v₀, integrator::EM_Kernel)
    f = forceUpdate(integrator.force, x₀)
    v_t=Vector{typeof(v₀)}()
    push!(v_t, deepcopy(v₀))
    noise=init_randn_correlated(integrator.σ_corr)
    return KernelEMState(deepcopy(x₀), deepcopy(v₀), f, v_t)
end



function UpdateState!(state::KernelEMState, integrator::EM_Kernel; kwargs...)
    mem_int= sum(integrator.kernel[:,:,i].*state.v_t[i] for i in 1:length(state.v_t); init=zeros(integrator.dim))
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)
    @. state.x = state.x + integrator.Δt * state.v
    apply_space!(integrator.bc,state.x,state.v)
    state.v = state.v .+integrator.Δt / integrator.M * state.f .+ integrator.Δt*mem_ker .+ randn_correlated(state,integrator)

    if length(state.v_t) == (size(kernel,3)) # In that case, we remove thing from start
        pop!(state.v_t)
    end
    pushfirst!(state.v_t, state.v)

    return nostop
end
