using LangevinDiffEq,DiffEqDevTools
using Plots

f1_harmonic(v,u,p,t) = -u
f2_harmonic(v,u,p,t) = v
γ = 1.0
g(u,p,t) = γ*exp(-t)

u0 = 0
v0 = 1

ff_harmonic = DynamicalSDEFunction(f1_harmonic,f2_harmonic,g)
prob1 = DynamicalSDEProblem(ff_harmonic,g,v0,u0,(0.0,5.0))

dts = (1/2) .^ (8:-1:4)

# Can't use NoiseGrid as noise is not generated with the correct size in convergence.jl. We require noise with shape of v.
sim1  = analyticless_test_convergence(dts,prob1,GJ(),(1/2)^10;trajectories=Int(2e2),use_noise_grid=false)
plot(sim1)
