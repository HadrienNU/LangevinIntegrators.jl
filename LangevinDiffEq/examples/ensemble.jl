using LangevinDiffEq
using DifferentialEquations.EnsembleAnalysis
using Plots

f1_harmonic(v,u,p,t) = -u
f2_harmonic(v,u,p,t) = v
γ = 1
g(u,p,t) = γ*ones(1)

u0 = zeros(1)
v0 = ones(1)

ff_harmonic = DynamicalSDEFunction(f1_harmonic,f2_harmonic,g)
prob = DynamicalSDEProblem(ff_harmonic,g,v0,u0,(0.0,10.0))
ensembleprob = EnsembleProblem(prob)

sol = solve(ensembleprob, GJ(), EnsembleThreads(), trajectories = 1000; dt=1/10)

summ = EnsembleSummary(sol, 0:0.1:10)
plot(summ, labels = "Middle 95%")
summ = EnsembleSummary(sol, 0:0.1:10; quantiles = [0.25, 0.75])
plot!(summ, labels = "Middle 50%", legend = true)
