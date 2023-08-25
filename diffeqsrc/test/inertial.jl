# En fait on va commencer par écrire les tests

# Tests existants pour BAOAB
using StochasticDiffEq, DiffEqNoiseProcess, Test, DiffEqDevTools, Random
Random.seed!(1)

f1_harmonic(v,u,p,t) = -u
f2_harmonic(v,u,p,t) = v
γ = 1
g(u,p,t) = γ # Ca devrait être gamma ? On doit poiuvoir dans un second implémenter une version de l'aglo qui a gamma constant


@testset "Scalar u" begin
    u0 = 0
    v0 = 1

    ff_harmonic = DynamicalSDEFunction(f1_harmonic,f2_harmonic,g)
    prob1 = DynamicalSDEProblem(ff_harmonic,ff_harmonic.g,v0,u0,(0.0,5.0))

    dts = (1/2) .^ (8:-1:4)

    # Can't use NoiseGrid as noise is not generated with the correct size in convergence.jl. We require noise with shape of v.
    sim1  = analyticless_test_convergence(dts,prob1,GJ(),(1/2)^10;trajectories=Int(2e2),use_noise_grid=false)
    display(sim1.𝒪est)
    @test abs(sim1.𝒪est[:weak_final]-1) < 0.3
end

@testset "Vector u" begin

    u0 = zeros(2)
    v0 = ones(2)

    f1_harmonic_iip(dv,v,u,p,t) = dv .= f1_harmonic(v,u,p,t)
    f2_harmonic_iip(du,v,u,p,t) = du .= f2_harmonic(v,u,p,t)
    g_iip(du,u,p,t) = du .= g(u,p,t)

    ff_harmonic = DynamicalSDEFunction(f1_harmonic,f2_harmonic,g)
    prob1 = DynamicalSDEProblem(ff_harmonic,g,v0,u0,(0.0,5.0))
    sol1 = solve(prob1,BAOAB(gamma=γ);dt=1/10,save_noise=true)

    prob2 = DynamicalSDEProblem(f1_harmonic_iip,f2_harmonic_iip,g_iip,v0,u0,(0.0,5.0); noise=NoiseWrapper(sol1.W))
    sol2 = solve(prob2,BAOAB(gamma=γ);dt=1/10)

    @test sol1[:] ≈ sol2[:]

    dts = (1/2) .^ (8:-1:4)

    # Can't use NoiseGrid as noise is not generated with the correct size in convergence.jl. We require noise with shape of v.
    sim1  = analyticless_test_convergence(dts,prob1,GJ(),(1/2)^10;trajectories=Int(1e2),use_noise_grid=false)
    @test abs(sim1.𝒪est[:weak_final]-1) < 0.3
end
