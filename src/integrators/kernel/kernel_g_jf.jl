struct GJF_Kernel{FP<:AbstractForce,TF<:AbstractFloat,TM} <: KernelIntegrator
    force::FP
    β::TF
    kernel::Array{TF}
    σ_corr::Array{TF}
    M::TM
    Δt::TF
    c₂::TF
    sc₁::TF
    d₁::TF
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
    a =  ker_mat[1,1,1] * Δt
    c₂ = (1 - 0.5* a) / (1 + 0.5 * a)
    sc₁ = sqrt((1+c₂)/2)
    d₁  = sqrt((1-c₂)/a)
    if dim==1
        noise_fdt=sqrt(M * Δt / β) * real.(ifft(sqrt.(fft(ker_mat,1)),1)) #Multiply value at 0 by 2?
    else
        noise_fdt=sqrt(M * Δt / β) * real.(ifft(sqrt.(fft(ker_mat,1)),1))# note quand Kernel est une matrix il faut faire le cholesky
    end
    return GJF_Kernel(force, β, ker_mat, noise_fdt, M, Δt, c₂, sc₁, d₁, dim, bc)
end



function UpdateState!(state::MemoryKernelState, integrator::GJF_Kernel; kwargs...)

    randn_correlated(state.ξ, integrator.σ_corr)
    memory_integral(state.memory, integrator)

    state.v_mid = integrator.sc₁ * state.v + 0.5* integrator.d₁ *integrator.Δt*state.f/integrator.M + 0.5*integrator.d₁*(state.ξ.int_val-state.memory.int_val)

    state.x = state.x .+ integrator.d₁*integrator.Δt * state.v_mid

    apply_space!(integrator.bc,state.x,state.v)
    nostop = forceUpdate!(integrator.force, state.f, state.x; kwargs...)

    state.v = (integrator.c₂ * state.v_mid +  0.5*integrator.d₁* integrator.Δt *state.f/integrator.M  + 0.5*integrator.d₁*(state.ξ.int_val-state.memory.int_val))/integrator.sc₁

    store_new_value(state.memory, state.v_mid)

    return nostop
end
