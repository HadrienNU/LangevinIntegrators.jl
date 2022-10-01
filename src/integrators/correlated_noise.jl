function init_randn_correlated(σ_corr::Array{TF}) where {TF<:AbstractFloat}
    return Vector{Vector{TF}}([randn(size(σ_corr,1)) for i in 1:size(σ_corr,3)])
end

function randn_correlated(state::AbstractMemoryKernelState, integrator::KernelIntegrator)
    # Compute the next correlated gaussian noise
    # Use fourrier transform to decompose the noise correlation when K(t-t')
    #Add and remove random number
    popfirst!(state.noise_n)
    push!(state.noise_n, randn(integrator.dim))
    return sum(integrator.σ_corr[:,:,i].*state.noise_n[i] for i in 1:size(integrator.σ_corr,1); init=zeros(integrator.dim))
end


function randn_correlated_karhunen_loeve(state::AbstractMemoryKernelState, integrator::KernelIntegrator)
    # Compute the next correlated gaussian noise
    # Use Karhunen–Loève théorem to decompose the noise correlation when K(t,t')
end
