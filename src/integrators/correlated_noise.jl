function init_randn_correlated(vector_size)
    noise=Queue{typeof(integrator.σ_corr[1])}()
    for n=1:vector_size
        push!(noise, randn(state.dim))
    end
    return noise
end

function randn_correlated(state::AbstractMemoryKernelState, integrator::KernelIntegrator)
    # Compute the next correlated gaussian noise
    # Use fourrier transform to decompose the noise correlation when K(t-t')
    #Add and remove random number
    popfirst!(state.noise_n)
    push!(state.noise_n, randn(state.dim))
    return sum(integrator.σ_corr.*state.noise_n)
end


function randn_correlated_karhunen_loeve(state::AbstractMemoryKernelState, integrator::KernelIntegrator)
    # Compute the next correlated gaussian noise
    # Use Karhunen–Loève théorem to decompose the noise correlation when K(t,t')
end


function corr_mem(v_t,kernel)
    #Compute the convolution of v_t and kernel
end
