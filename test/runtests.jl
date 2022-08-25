using LangevinIntegrators
using Test


@testset "LangevinIntegrators.jl" begin

    @testset "Read Ini conf" begin
        # Test from_config with various config file
        params, init_conds_args = read_conf("test_params.ini")

        @test params.n_steps == 5000
        @test params.n_trajs == 1

        integrator = read_integrator_conf("test_integrator_em.ini")

        cond_arr = initialize_initcond(integrator, init_conds_args)

        @test length(cond_arr) == 1
        @test cond_arr[1] isa LangevinIntegrators.Gaussian_InitCond
        @test cond_arr[1].mean == [0.0]


        integrator = read_integrator_conf("test_integrator_baoab.ini")

        cond_arr = initialize_initcond(integrator, init_conds_args)

        @test length(cond_arr) == 2
        @test cond_arr[2] isa LangevinIntegrators.Gaussian_InitCond
        @test cond_arr[2].mean == [0.0]
        @test cond_arr[2].std == 1.0

        # Test that the force from the config is correct
        @test force_eval(integrator.force, [0.0])[1] ≈ 0.0
        @test force_eval(integrator.force, [1.0])[1] ≈ -2.0


        params, init_conds_args = read_conf("test_hidden.ini")
        integrator = read_integrator_conf("test_hidden.ini")

        cond_arr = initialize_initcond(integrator, init_conds_args)

        @test length(cond_arr) == 3
        @test cond_arr[3] isa LangevinIntegrators.Gaussian_InitCond
        @test cond_arr[3].mean ≈ [-0.157145581734742, -0.0034967990480131105]
        @test cond_arr[3].std == 1.0


    end

    @testset "NPZread" begin
        # Test from npz

        integrator = read_integrator_hidden_npz("test_coeffs.npz")

        @test integrator.Δt ≈ 0.001

        @test force_eval(integrator.force, [0.5])[1] ≈ -0.07372446770411359
        @test force_eval(integrator.force, [1.0])[1] ≈ -9.0294072


        int_free_energy = read_integrator_hidden_npz("test_free_energy_force.npz")

        @test force_eval(int_free_energy.force, [1.0])[1] ≈ 49.65808451
        @test force_eval(int_free_energy.force, [2.0])[1] ≈ -0.3815412

        int_linear = read_integrator_hidden_npz("test_linear_force.npz")

        @test int_linear.dim_tot == 2
        @test size(int_linear.S) == (2, 2)
        @test int_linear.friction_hv[1, 1] ≈ -0.024106487767047175
        @test int_linear.friction_hh[1, 1] ≈ 0.0147817703

        @test force_eval(int_linear.force, [0.0])[1] ≈ 0.0
        @test force_eval(int_linear.force, [1.0])[1] ≈ -1.00027734

    end

    @testset "Forces" begin
        # From potential
        force = ForceFromPotential("Harmonic")
        @test force_eval(force, [0.0])[1] ≈ 0.0
        @test force_eval(force, [1.0])[1] ≈ -1.0

        force = ForceFromPotential("DoubleWell")
        @test force_eval(force, [0.0])[1] ≈ 0.0
        @test force_eval(force, [1.0])[1] ≈ 0.0

        # force=ForceFromPotential("Muller")
        # @test force_eval(force,[0.0 0.0]) ≈ 0.0
        # @test force_eval(force,[1.0 1.0])[1] ≈ 0.0

        #Forces from ApproxFun
        force = ForceFromBasis("Taylor", [0.0 -2.0])
        @test force_eval(force, [0.0])[1] ≈ 0.0
        @test force_eval(force, [1.0])[1] ≈ -2.0

        #Forces from BSplinesKit
        #[0.163005, 0.163005, 0.163005, 0.163005, 0.622277, 1.081549, 1.081549, 1.081549, 1.081549]
        force = ForceFromSplines(
            3,
            [0.163005, 0.163005, 0.163005, 0.163005, 0.622277, 1.081549, 1.081549, 1.081549, 1.081549],
            [8.80178094 -0.25194201 1.65292369 -8.13424986 -10.49673538],
        )

        @test force_eval(force, [0.5])[1] ≈ -0.07372447
        @test force_eval(force, [1.0])[1] ≈ -9.0294072
        # @test force_eval(force,[2.0])[1] ≈ 49.3164884 # BSplineKit does not extrapolate TODO

        #Forces from Scipy.interpolate
        force = ForceFromScipySplines(
            3,
            [0.163005, 0.163005, 0.163005, 0.163005, 0.622277, 1.081549, 1.081549, 1.081549, 1.081549],
            [8.80178094 -0.25194201 1.65292369 -8.13424986 -10.49673538],
        )

        @test force_eval(force, [0.5])[1] ≈ -0.07372447
        @test force_eval(force, [1.0])[1] ≈ -9.0294072
        @test force_eval(force, [2.0])[1] ≈ 49.3164884
    end

    @testset "Initial conditions" begin # Voir si on peut faire une matrice de test

        params_full, init_conds_args = read_conf("test_hidden.ini")

        #Tester les différents moyens d'avoir des conditions initiales
        force = ForceFromPotential("Harmonic")
        #Overdamped
        integrator = EM(force, 1.0, 1e-3)
        cond_arr = initialize_initcond(integrator, init_conds_args)
        state = InitState(integrator, cond_arr)
        @test state.x ≈ [1.0]
        @test state isa LangevinIntegrators.AbstractOverdampedState

        #Underdamped
        integrator = Verlet(force, 1.0, 1e-3)
        @test state.x ≈ [1.0]
        cond_arr = initialize_initcond(integrator, init_conds_args)
        state = InitState(integrator, cond_arr)
        @test state isa LangevinIntegrators.AbstractInertialState

        #Hidden
        integrator = EM_Hidden(force, [[1.0, 1.0] [-1.0, 2.0]], [[1.0, 0.0] [0.0, 1.0]], 1e-3, 1)
        cond_arr = initialize_initcond(integrator, init_conds_args)
        state = InitState(integrator, cond_arr)
        @test state.x ≈ [1.0]
        @test state isa LangevinIntegrators.AbstractMemoryHiddenState
    end


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

    @testset "run_trajectory" begin # And then a matrix over the integrator

        init_conds_args = Dict("position" => Dict("type" => "Cste"), "velocity" => Dict("type" => "Cste"))
        force = ForceFromPotential("Harmonic")
        params = TrajsParams(; n_steps = 10^3)

        integrator = EM(force, 1.0, 0.001)
        init_conds = initialize_initcond(integrator, init_conds_args)

        state = InitState(integrator, init_conds)
        time, stop = run_trajectory!(state, integrator; params = params)

        # We should do statiscal test
        @test time ≈ params.n_steps * integrator.Δt
        @test stop == false

    end

    @testset "run_multiple_trajectories" begin
        params = TrajsParams(; n_steps = 10^3, n_trajs = 10)
        init_conds_args = Dict("position" => Dict("type" => "Cste"), "velocity" => Dict("type" => "Cste"))

        force = ForceFromPotential("Harmonic")
        integrator = EM(force, 1.0, 0.001)

        run_trajectories(integrator; params = params, init_conds_args = init_conds_args)
    end

    include("plumed_tests.jl")

end
