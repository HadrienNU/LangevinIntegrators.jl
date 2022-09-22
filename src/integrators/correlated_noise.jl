function randn_correlated(state::AbstractMemoryKernelState, integrator::KernelIntegrator)
    # Compute the next correlated gaussian noise
    # Use fourrier transform to decompose the noise correlation when K(t-t')
end


function randn_correlated_karhunen_loeve(state::AbstractMemoryKernelState, integrator::KernelIntegrator)
    # Compute the next correlated gaussian noise
    # Use Karhunen–Loève théorem to decompose the noise correlation when K(t,t')
end


function corr_mem(v_t,kernel)
    #Compute the convolution of v_t and kernel
end
