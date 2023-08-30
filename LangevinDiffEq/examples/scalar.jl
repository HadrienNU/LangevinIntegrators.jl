using LangevinDiffEq
using Plots

f1_harmonic(v,u,p,t) = -u
f2_harmonic(v,u,p,t) = v
γ = 1.0
g(u,p,t) = γ*exp(-t)

u0 = 0
v0 = 1

ff_harmonic = DynamicalSDEFunction(f1_harmonic,f2_harmonic,g);
prob = DynamicalSDEProblem(ff_harmonic,g,v0,u0,(0.0,5.0)); # , noise_rate_prototype=γ
sol = solve(prob,GJ();dt=1/10);


plot(sol)
