using LangevinDiffEq,DiffEqNoiseProcess
using Plots

γ = 1
u0 = zeros(1)
v0 = ones(1)



f1_harmonic(v,u,p,t) = -u
f2_harmonic(v,u,p,t) = v
g(u,p,t) = γ*exp(-t).*ones(1)

ff_harmonic = DynamicalSDEFunction(f1_harmonic,f2_harmonic,g);
prob1 = DynamicalSDEProblem(ff_harmonic,g,v0,u0,(0.0,5.0)); # , noise_rate_prototype=γ
sol1 = solve(prob1,GJ();dt=1/10,save_noise=true);


println("In place case")
f1_harmonic_iip(dv,v,u,p,t) = dv .= f1_harmonic(v,u,p,t)
f2_harmonic_iip(du,v,u,p,t) = du .= f2_harmonic(v,u,p,t)
g_iip(du,u,p,t) = du .=  g(u,p,t)

ff_harmonic = DynamicalSDEFunction(f1_harmonic_iip,f2_harmonic_iip,g_iip)
prob = DynamicalSDEProblem(ff_harmonic,g_iip,v0,u0,(0.0,5.0); noise=NoiseWrapper(sol1.W))
sol = solve(prob,GJ();dt=1/10);

plot(sol1)
plot!(sol)
