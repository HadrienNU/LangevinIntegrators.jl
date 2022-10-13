struct BBK_Kernel{FP<:AbstractForce,TF<:AbstractFloat,TM} <: KernelIntegrator
    force::FP
    β::TF
    kernel::Array{TF}
    σ_corr::Array{TF}
    M::TM
    Δt::TF
    invK::Union{TF,Matrix{TF}}
    dim::Int64
    bc::Union{AbstractSpace,Nothing}
end

"""
    BBK_Kernel(force, β, γ, M, Δt)

Set up the BBK_Kernel integrator for generalized Langevin.

Adapted from Iterative Reconstruction of Memory Kernels Gerhard Jung,*,†,‡ Martin Hanke,*,§ and Friederike Schmid*,†

### Fields

* force   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
function BBK_Kernel(force::FP, β::TF, kernel::Array{TF}, M::TM, Δt::TF, dim::Int64=1, bc::Union{AbstractSpace,Nothing}=nothing) where {FP<:AbstractForce,TF<:AbstractFloat,TM}
    if  size(kernel,2) != dim || size(kernel,2) != size(kernel,3)
        throw(ArgumentError("Mismatch of dimension bewteen kernel and space dimension"))
    end
    ker_mat=reshape(kernel,:,dim,dim)
    invK = inv(1 .+ 0.5 * Δt * ker_mat[1,:,:])
    if dim==1
        noise_fdt=sqrt(Δt / β) * real.(ifft(sqrt.(fft(ker_mat,1)),1))
    else
        noise_fdt=sqrt(Δt / β) * real.(ifft(sqrt.(fft(ker_mat,1)),1))# note quand Kernel est une matrix il faut faire le cholesky
    end
    return BBK_Kernel(force, β, ker_mat, noise_fdt, M, Δt, invK, dim, bc)
end

mutable struct BBKKernelState{TF<:AbstractFloat} <: AbstractMemoryKernelState
    x::Vector{TF}
    v::Vector{TF}
    v_mid::Vector{TF}
    f::Vector{TF}
    diss_f::Vector{TF} # Dissipative force
    v_t::Vector{Vector{TF}} # Trajectory of v to compute the kernel
    noise_n::Vector{Vector{TF}} # Utiliser un Circular buffer à la place?

    function BBKKernelState(x₀::Vector{TF}, v₀::Vector{TF}, f::Vector{TF}, v_t, noise::Vector{Vector{TF}}) where {TF<:AbstractFloat}
        return new{TF}(x₀, v₀, similar(v₀), f, similar(f) , v_t , noise) # TODO: import v_t as a vector and fill the Vector by zeros and given v_t
    end
end

function InitState!(x₀, v₀, integrator::BBK_Kernel)
    f = forceUpdate(integrator.force, x₀)
    v_t=Vector{typeof(v₀)}()
    push!(v_t, v₀)
    noise=init_randn_correlated(integrator.σ_corr)
    return BBKKernelState(x₀, v₀, f, v_t, noise)
end

function UpdateState!(state::BBKKernelState, integrator::BBK_Kernel; kwargs...)
    state.v_mid = state.v + 0.5 * integrator.Δt / integrator.M * state.f .-integrator.Δt .*( state.diss_f+ 0.5*integrator.kernel[1,:,:]*state.v)
    @. state.x = state.x + integrator.Δt * state.v_mid
    apply_space!(integrator.bc,state.x,state.v)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)

    state.diss_f = zeros(integrator.dim)
    for k in 1:integrator.dim, l in 1:integrator.dim
        for i in 2:size(state.v_t,1)
            @inbounds state.diss_f[k] += integrator.kernel[i,k,l]*state.v_t[i][l]
        end
    end
    state.diss_f +=randn_correlated(state, integrator)
    # state.diss_f = sum(integrator.kernel[:,:,i]*state.v_t[i] for i in 2:length(state.v_t); init=zeros(integrator.dim)) + randn_correlated(state, integrator) # Pour ne calculer l'intégrale qu'une fois, on la stocke puisqu'elle resservira au prochain pas de temps
    state.v = integrator.invK * (state.v_mid + 0.5 * integrator.Δt / integrator.M * state.f + integrator.Δt*state.diss_f)

    if length(state.v_t) == (size(integrator.kernel,1)) # In that case, we remove thing from start
        pop!(state.v_t)
    end
    pushfirst!(state.v_t, state.v)

    return nostop
end
