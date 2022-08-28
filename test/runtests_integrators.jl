
@testset "Boundary conditions" begin # Voir si on peut faire une matrice de test

end

@testset "integrators_overdamped" begin
    force = ForceFromPotential("Harmonic")

    integrator = EM(force, 1.0, 1e-3)
    state = InitState!([0.0], integrator)
    UpdateState!(state, integrator)

end


@testset "integrators_inertial" begin
    force = ForceFromPotential("Harmonic")

    integrator = BBK(force, 1.0, 1.0, 1.0, 1e-3)
    state = InitState!([0.0], [0.0], integrator)
    UpdateState!(state, integrator)

    integrator = GJF(force, 1.0, 1.0, 1.0, 1e-3)
    state = InitState!([0.0], [0.0], integrator)
    UpdateState!(state, integrator)

    integrator = ABOBA(force, 1.0, 1.0, 1.0, 1e-3)
    state = InitState!([0.0], [0.0], integrator)
    UpdateState!(state, integrator)

    integrator = BAOAB(force, 1.0, 1.0, 1.0, 1e-3)
    state = InitState!([0.0], [0.0], integrator)
    UpdateState!(state, integrator)

    integrator = Verlet(force, 1.0, 1e-3)
    state = InitState!([0.0], [0.0], integrator)
    UpdateState!(state, integrator)
    #Comme Verlet est déterministe on peut faire des test sur la valeur
    # @test state.x ≈
end

@testset "integrators_hidden" begin
    force = ForceFromPotential("Harmonic")

    integrator = EM_Hidden(force, [[1.0, 1.0] [-1.0, 2.0]], [[1.0, 0.0] [0.0, 1.0]], 1e-3, 1)
    state = InitState!([0.0], [0.0], [0.0], integrator)
    UpdateState!(state, integrator)

    integrator = ABOBA_Hidden(force, [[1.0, 1.0] [-1.0, 2.0]], [[1.0, 0.0] [0.0, 1.0]], 1e-3, 1)
    state = InitState!([0.0], [0.0], [0.0], integrator)
    UpdateState!(state, integrator)

end

@testset "run_trajectory overdamped $int_struct" for int_struct in [EM] # And then a matrix over the integrator

    init_conds_args = Dict("position" => Dict("type" => "Cste"), "velocity" => Dict("type" => "Cste"))
    force = ForceFromPotential("Harmonic")
    params = TrajsParams(; n_steps = 10^3)

    integrator = int_struct(force, 1.0, 0.001)
    init_conds = initialize_initcond(integrator, init_conds_args)

    state = InitState(integrator, init_conds)
    time, stop = run_trajectory!(state, integrator; params = params)

    # We should do statiscal test
    @test time ≈ params.n_steps * integrator.Δt
    @test stop == false

end

@testset "run_trajectory inertial $int_struct" for int_struct in [ABOBA,BAOAB,BBK,GJF] # And then a matrix over the integrator

    init_conds_args = Dict("position" => Dict("type" => "Cste"), "velocity" => Dict("type" => "Cste"))
    force = ForceFromPotential("Harmonic")
    params = TrajsParams(; n_steps = 10^3)

    integrator = int_struct(force, 1.0,0.1,2.0, 0.001)
    init_conds = initialize_initcond(integrator, init_conds_args)

    state = InitState(integrator, init_conds)
    time, stop = run_trajectory!(state, integrator; params = params)

    # We should do statiscal test
    @test time ≈ params.n_steps * integrator.Δt
    @test stop == false

end
