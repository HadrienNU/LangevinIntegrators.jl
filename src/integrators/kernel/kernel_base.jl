mutable struct CircularVector{TF<:AbstractFloat} # TODO replace the vector by a CircularBuffer?
    values::Vector{Vector{TF}}
    int_val::Vector{TF}
    curr_ind::Int64
    len::Int64
    dim::Int64
    function CircularVector(values::Vector{Vector{TF}}) where {TF<:AbstractFloat}
        if length(values) <1
            throw(ArgumentError("CircularVector cannot be not empty"))
        end
        return new{TF}(values,similar(values[1]),1,length(values),length(values[1]))
    end
end

function randomCircularVector(σ_corr::Array{TF}) where {TF<:AbstractFloat}
    return CircularVector([randn(size(σ_corr,2)) for i in 1:size(σ_corr,1)])
end

function zerosCircularVector(len::Int64,dim::Int64, TF)
    return CircularVector([zeros(TF,dim) for i in 1:len])
end

function Base.length(vec::CircularVector)
    return vec.len
end

function Base.getindex(vec::CircularVector,ind::Int)
    return vec.values[1+(ind+vec.curr_ind-1)%length(vec)]
end


function store_new_value(vec::CircularVector, val::Vector{TF}) where {TF<:AbstractFloat}
    vec.curr_ind -= 1
    if vec.curr_ind == 0
        vec.curr_ind = length(vec)
    end
    vec.values[vec.curr_ind] = val
end


function randn_correlated(noise_vec::CircularVector, σ_corr::Array{TF}) where {TF<:AbstractFloat}
    # Compute the next correlated gaussian noise
    store_new_value(noise_vec,randn(noise_vec.dim))
    noise_vec.int_val[:] .= 0. #zeros(noise_vec.dim)
    for l in 1:noise_vec.dim, k in 1:noise_vec.dim
        for i in 1:size(σ_corr,1)
            @inbounds noise_vec.int_val[k] += σ_corr[i,k,l]*noise_vec[i][l]
        end
    end
    return noise_vec.int_val
end

function randn_correlated_karhunen_loeve(state::AbstractMemoryKernelState, integrator::KernelIntegrator)
    # Compute the next correlated gaussian noise
    # Use Karhunen–Loève théorem to decompose the noise correlation when K(t,t')
end


function memory_integral(vt::CircularVector, integrator::MKI) where {MKI <: KernelIntegrator} # Keep integrator as argument as it allow to specialize the function if needed for various integrator
    vt.int_val[:] .= 0. #zeros(vt.dim)
    for l in 1:vt.dim, k in 1:vt.dim
        for i in 2:size(integrator.kernel,1)
            @inbounds vt.int_val[k] += integrator.kernel[i,k,l]*vt[i][l]
        end
    end
    vt.int_val*=integrator.Δt
    return vt.int_val
end

function memory_integral(vt::CircularVector, kernel::Array{TF}, Δt::TF) where {TF<:AbstractFloat}
    vt.int_val[:] .= 0. # zeros(vt.dim)
    for l in 1:vt.dim, k in 1:vt.dim
        for i in 2:size(kernel,1)
            @inbounds vt.int_val[k] += kernel[i,k,l]*vt[i][l]* Δt
        end
    end
    return vt.int_val
end


mutable struct MemoryKernelState{TF<:AbstractFloat} <: AbstractMemoryKernelState
    x::Vector{TF}
    v::Vector{TF}
    v_mid::Vector{TF}
    f::Vector{TF}
    memory::CircularVector{TF} # Trajectory of v to compute the kernel, both array are given by the size of the kernel
    ξ::CircularVector{TF}
    function MemoryKernelState(x₀::Vector{TF}, v₀::Vector{TF}, f::Vector{TF}, σ_corr::Array{TF}) where {TF<:AbstractFloat}
        ξ = randomCircularVector(σ_corr)
        randn_correlated(ξ,σ_corr)
        memory = zerosCircularVector(size(σ_corr,1),length(v₀), TF)
        store_new_value(memory, v₀)
        return new{TF}(x₀, v₀, similar(v₀), f, memory, ξ)
    end
end

mutable struct MemoryKernelStateTwoNoise{TF<:AbstractFloat} <: AbstractMemoryKernelState
    x::Vector{TF}
    v::Vector{TF}
    v_mid::Vector{TF}
    f::Vector{TF}
    memory::CircularVector{TF} # Trajectory of v to compute the kernel, both array are given by the size of the kernel
    ξ::CircularVector{TF}
    ξ₂::CircularVector{TF}
    function MemoryKernelStateTwoNoise(x₀::Vector{TF}, v₀::Vector{TF}, f::Vector{TF}, σ_corr::Array{TF}) where {TF<:AbstractFloat}
        ξ = randomCircularVector(σ_corr)
        randn_correlated(ξ,σ_corr)
        ξ₂=randomCircularVector(σ_corr)
        randn_correlated(ξ₂,σ_corr)
        memory = zerosCircularVector(size(σ_corr,1),length(v₀), TF)
        store_new_value(memory,v₀)
        return new{TF}(x₀, v₀, similar(v₀), f, memory, ξ, ξ₂)
    end
end

function InitState!(x₀, v₀, integrator:: MKI ) where {MKI <: KernelIntegrator}
    if integrator.dim != length(x₀)
        throw(ArgumentError("Mismatch of dimension in state initialization"))
    end
    f = forceUpdate(integrator.force, x₀)
    return MemoryKernelState(x₀, v₀, f, integrator.σ_corr)
end
