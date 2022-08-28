using LangevinIntegrators
using Test
using DelimitedFiles

@testset "LangevinIntegrators.jl" begin

    include("runtests_init.jl")

    include("runtests_forces.jl")

    include("runtests_integrator.jl")


    @testset "run_multiple_trajectories" begin
        params = TrajsParams(; n_steps = 10^3, n_trajs = 10)
        init_conds_args = Dict("position" => Dict("type" => "Cste"), "velocity" => Dict("type" => "Cste"))

        force = ForceFromPotential("Harmonic")
        integrator = EM(force, 1.0, 0.001)

        run_trajectories(integrator; params = params, init_conds_args = init_conds_args)
    end


    @testset "observers" begin # And then a matrix over the integrator

        init_conds_args = Dict("position" => Dict("type" => "Cste"), "velocity" => Dict("type" => "Cste"))
        force = ForceFromPotential("Harmonic")
        params = TrajsParams(; n_steps = 10^3)

        integrator = EM(force, 1.0, 0.001)
        init_conds = initialize_initcond(integrator, init_conds_args)

        state = InitState(integrator, init_conds)
        time, stop = run_trajectory!(state, integrator; params = params)

        # @test isfile("trajectory.dat")
        # rm("trajectory.dat")
    end

    include("plumed_tests.jl")

end
