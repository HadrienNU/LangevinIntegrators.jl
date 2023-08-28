using LangevinDiffEq

γ = 1
u0 = zeros(2)
v0 = ones(2)

f1_harmonic_iip(dv,v,u,p,t) = dv .= -u
f2_harmonic_iip(du,v,u,p,t) = du .= v
g_iip(du,u,p,t) = du .=  γ

ff_harmonic = DynamicalSDEFunction(f1_harmonic_iip,f2_harmonic_iip,g_iip)
prob1 = DynamicalSDEProblem(ff_harmonic,g_iip,v0,u0,(0.0,5.0))
sol1 = solve(prob1,GJ();dt=1/10)
