using LangevinIntegrators
using Test
using DelimitedFiles

@testset "LangevinIntegrators.jl" begin

    include("runtests_init.jl")

    include("runtests_forces.jl")

    include("runtests_integrators.jl")


    @testset "run_multiple_trajectories" begin
        params = TrajsParams(; n_steps = 10^3, n_trajs = 5, n_save_iters=49, save_filename_pattern = "trajectory_*.dat") # Ajouter des paramètres pour sauvegarder les trajectoires
        init_conds_args = Dict("position" => Dict("type" => "Cste"), "velocity" => Dict("type" => "Cste"))

        force = ForceFromPotential("Harmonic")
        integrator = EM(force, 1.0, 0.001)

        trajs=run_trajectories(integrator; params = params, init_conds_args = init_conds_args, buffer_size = 7)
        list_files=filter(s -> occursin(r"trajectory_.*.dat",s),readdir())

        @test length(list_files) == params.n_trajs
        for n in 1:length(list_files)
            data=readdlm(list_files[n])
            @test data[end,1] ≈ trajs[n].time[trajs[n].save_index]
            @test data[end,2:end] ≈ trajs[n].xt[trajs[n].save_index,:]

            @test trajs[n].buffer_size == 7
            @test trajs[n].save_index == 6
            rm(list_files[n]) # clean behind test
        end

        fpt, reached=run_fpt(integrator; params = params, init_conds_args = init_conds_args)

        @test length(fpt) == params.n_trajs
        @test fpt[end] > 0.0
        @test sum(reached) == 0 # None of them should have reached, since thre is no stop conditions here

        @test sum(fpt) / params.n_trajs ≈ params.n_steps * integrator.Δt

    end


    # @testset "observers" begin # And then a matrix over the integrator
    #
    #     init_conds_args = Dict("position" => Dict("type" => "Cste"), "velocity" => Dict("type" => "Cste"))
    #     force = ForceFromPotential("Harmonic")
    #     params = TrajsParams(; n_steps = 10^3)
    #
    #     integrator = EM(force, 1.0, 0.001)
    #     init_conds = initialize_initcond(integrator, init_conds_args)
    #
    #     state = InitState(integrator, init_conds)
    #     time, stop = run_trajectory!(state, integrator; params = params)
    #
    #     # @test isfile("trajectory.dat")
    #     # rm("trajectory.dat")
    # end

    include("plumed_tests.jl")

end
