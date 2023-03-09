struct BBK_Kernel{FP<:AbstractForce,TF<:AbstractFloat,TM} <: KernelIntegrator
    force::FP
    β::TF
    kernel::Array{TF}
    σ_corr::Array{TF}
    M::TM
    Δt::TF
    c₀::Union{TF,Matrix{TF}}
    c₁::Union{TF,Matrix{TF}}
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
    c₀ = 1 .- 0.5 * Δt * ker_mat[1,:,:]
    return BBK_Kernel(force, β, ker_mat, noise_fdt, M, Δt, c₀, invK, dim, bc)
end


function UpdateState!(state::MemoryKernelState, integrator::BBK_Kernel; kwargs...)

    state.v_mid = integrator.c₀ * state.v .+ 0.5 * integrator.Δt / integrator.M * state.f  .- 0.5* integrator.Δt*state.memory.int_val .+ 0.5*state.ξ.int_val
    @. state.x = state.x + integrator.Δt * state.v_mid

    apply_space!(integrator.bc,state.x,state.v)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)

    randn_correlated(state.ξ, integrator.σ_corr)
    memory_integral(state.memory, integrator)

    state.v = integrator.c₁ * (state.v_mid .+ 0.5 * integrator.Δt / integrator.M * state.f  - 0.5* integrator.Δt*state.memory.int_val + 0.5*state.ξ.int_val)

    store_new_value(state.memory,state.v)

    return nostop
end
