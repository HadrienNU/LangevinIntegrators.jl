struct GJF_Kernel{FP<:AbstractForce,TF<:AbstractFloat,TM} <: KernelIntegrator
    force::FP
    β::TF
    kernel::Array{TF}
    σ_corr::Array{TF}
    M::TM
    Δt::TF
    a::Union{TF,Matrix{TF}}
    b::Union{TF,Matrix{TF}}
    dim::Int64
    bc::Union{AbstractSpace,Nothing}
end

"""
    GJF_Kernel(force, β, γ, M, Δt)

Set up the G-JF integrator for inertial Langevin.

Adapted from Iterative Reconstruction of Memory Kernels Gerhard Jung,*,†,‡ Martin Hanke,*,§ and Friederike Schmid*,†

### Fields

* force   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
function GJF_Kernel(force::FP, β::TF, kernel::Array{TF}, M::TM, Δt::TF, dim::Int64=1, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
    if  size(kernel,2) != dim || size(kernel,2) != size(kernel,3)
        throw(ArgumentError("Mismatch of dimension bewteen kernel and space dimension"))
    end
    ker_mat=reshape(kernel,:,dim,dim)

    a = (1 .- 0.5 * ker_mat[1,:,:] * Δt) *inv( (1 .+ 0.5 * ker_mat[1,:,:] * Δt))
    b = inv(1 .+ 0.5 * ker_mat[1,:,:] * Δt)
    if dim==1
        noise_fdt=sqrt(M * Δt / β) * real.(ifft(sqrt.(fft(ker_mat,1)),1)) #Multiply value at 0 by 2?
    else
        noise_fdt=sqrt(M * Δt / β) * real.(ifft(sqrt.(fft(ker_mat,1)),1))# note quand Kernel est une matrix il faut faire le cholesky
    end
    return GJF_Kernel(force, β, ker_mat, noise_fdt, M, Δt, a, b, dim, bc)
end

mutable struct GJFKernelState{TF<:AbstractFloat} <: AbstractMemoryKernelState
    x::Vector{TF}
    v::Vector{TF}
    f::Vector{TF}
    f_new::Vector{TF}
    ξ::Vector{TF}
    x_t::Vector{Vector{TF}} # Trajectory of x to compute the kernel, both array are given by the size of the kernel
    noise_n::Vector{Vector{TF}}
    function GJFKernelState(x₀::Vector{TF}, v₀::Vector{TF}, f::Vector{TF}, x_t, noise::Vector{Vector{TF}}) where {TF<:AbstractFloat}
        return new{TF}(x₀, v₀, f, copy(f), similar(f), x_t, noise)
    end
end

function InitState!(x₀, v₀, integrator::GJF_Kernel)
    f = forceUpdate(integrator.force, x₀)
    x_t=Vector{typeof(v₀)}()
    push!(x_t, x₀)
    noise=init_randn_correlated(integrator.σ_corr)
    return GJFKernelState(x₀, v₀, f, x_t, noise)
end



function UpdateState!(state::GJFKernelState, integrator::GJF_Kernel; kwargs...)

    state.ξ = randn_correlated(state,integrator)
    mem_int = zeros(integrator.dim)
    for l in 1:integrator.dim, k in 1:integrator.dim
        for i in 2:(length(state.x_t)-1)
            @inbounds mem_int[k] += integrator.kernel[i,k,l]*(state.x_t[i][l] - state.x_t[i-1][l])
        end
    end

    # mem_int = sum(integrator.kernel[:,:,i]*(state.x_t[i] - state.x_t[i-1]) for i in 2:(length(state.x_t)-1); init=zeros(integrator.dim))

    state.x = state.x .+ integrator.Δt .* integrator.b *state.v .+ 0.5 * (integrator.Δt^2 / integrator.M) .* integrator.b * state.f.- 0.5 * integrator.Δt * integrator.b * mem_int .+ 0.5 * (integrator.Δt/ integrator.M) .* integrator.b * state.ξ

    apply_space!(integrator.bc,state.x,state.v)
    nostop = forceUpdate!(integrator.force, state.f_new, state.x; kwargs...)

    state.v = integrator.a * state.v .+ 0.5 * integrator.Δt / integrator.M * (integrator.a * state.f .+ state.f_new).- integrator.b * mem_int .+ (1/ integrator.M) * integrator.b * state.ξ

    state.f = state.f_new

    if length(state.x_t) == (size(integrator.kernel,1)+1) # In that case, we remove thing from start
        pop!(state.x_t)
    end
    pushfirst!(state.x_t, state.x)

    return nostop
end
