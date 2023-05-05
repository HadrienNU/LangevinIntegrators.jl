
abstract type AbstractNoise end


struct GaussianNoise <: AbstractNoise
    dim::Int64
    function GaussianNoise(dim::Int64 = 1)
        return new(dim)
    end
end

function generate_noise!(ξ::Vector{TF},noise::GaussianNoise, x::Vector{TF},v::Vector{TF}) where {TF<:AbstractFloat}
    ξ .= randn(noise.dim)
end
#
# abstract type AbstractHistoryNoise end
#
# struct HistoryNoise where {TF<:AbstractFloat} <: AbstractHistoryNoise
#     dim::Int64
#     kernel::Array{TF}
#     σ_corr::Array{TF}
#     memory::CircularVector{TF} # Trajectory of v to compute the kernel, both array are given by the size of the kernel
#     ξ::CircularVector{TF}
# end
#
# function HistoryNoise(args)
#     body
# end
#
# function randn_correlated(
#     noise_vec::CircularVector,
#     σ_corr::Array{TF},
# ) where {TF<:AbstractFloat}
#     # Compute the next correlated gaussian noise
#     store_new_value(noise_vec, randn(noise_vec.dim))
#     noise_vec.int_val[:] .= 0.0 #zeros(noise_vec.dim)
#     for l = 1:noise_vec.dim, k = 1:noise_vec.dim
#         for i = 1:size(σ_corr, 1)
#             @inbounds noise_vec.int_val[k] += σ_corr[i, k, l] * noise_vec[i][l]
#         end
#     end
#     return noise_vec.int_val
# end
#
# function randn_correlated_karhunen_loeve(
#     state::AbstractMemoryKernelState,
#     integrator::KernelIntegrator,
# )
#     # Compute the next correlated gaussian noise
#     # Use Karhunen–Loève théorem to decompose the noise correlation when K(t,t')
# end
#
#
# function memory_integral(vt::CircularVector, integrator::MKI) where {MKI<:KernelIntegrator} # Keep integrator as argument as it allow to specialize the function if needed for various integrator
#     vt.int_val[:] .= 0.0 #zeros(vt.dim)
#     for l = 1:vt.dim, k = 1:vt.dim
#         for i = 2:size(integrator.kernel, 1)
#             @inbounds vt.int_val[k] += integrator.kernel[i, k, l] * vt[i][l]
#         end
#     end
#     vt.int_val *= integrator.Δt
#     return vt.int_val
# end
#
# function memory_integral(
#     vt::CircularVector,
#     kernel::Array{TF},
#     Δt::TF,
# ) where {TF<:AbstractFloat}
#     vt.int_val[:] .= 0.0 # zeros(vt.dim)
#     for l = 1:vt.dim, k = 1:vt.dim
#         for i = 2:size(kernel, 1)
#             @inbounds vt.int_val[k] += kernel[i, k, l] * vt[i][l] * Δt
#         end
#     end
#     return vt.int_val
# end
#
# function generate_noise!(ξ::Vector{TF},noise::HistoryNoise, x::::Vector{TF},v::::Vector{TF})
#     ξ.= randn_correlated(noise) .- memory_integral(noise)
# end
#
# abstract type AbstractHiddenNoise end
#
# struct HiddenNoise  <: AbstractHiddenoise
#     dim::Int64
#     kernel::Array{TF}
#     σ_corr::Array{TF}
#     memory::CircularVector{TF} # Trajectory of v to compute the kernel, both array are given by the size of the kernel
#     ξ::CircularVector{TF}
# end
