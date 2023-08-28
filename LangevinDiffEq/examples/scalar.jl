using LangevinDiffEq

f1_harmonic(v,u,p,t) = -u
f2_harmonic(v,u,p,t) = v
γ = 1.0
g(u,p,t) = γ

u0 = 0
v0 = 1

ff_harmonic = DynamicalSDEFunction(f1_harmonic,f2_harmonic,g)
prob1 = DynamicalSDEProblem(ff_harmonic,g,v0,u0,(0.0,5.0), noise_rate_prototype=γ)
sol1 = solve(prob1,GJ();dt=1/10)
