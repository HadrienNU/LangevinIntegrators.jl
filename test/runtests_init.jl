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