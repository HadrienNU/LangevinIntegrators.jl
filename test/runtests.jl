using LangevinIntegrators
using Test
using DelimitedFiles

@testset "LangevinIntegrators.jl" begin

    include("runtests_init.jl")

    include("runtests_forces.jl")

    include("runtests_integrator.jl")


    @testset "run_multiple_trajectories" begin
        params = TrajsParams(; n_steps = 10^3, n_trajs = 10, n_save_iters=10^2, save_filename_pattern = "trajectory_*.dat") # Ajouter des paramètres pour sauvegarder les trajectoires
        init_conds_args = Dict("position" => Dict("type" => "Cste"), "velocity" => Dict("type" => "Cste"))

        force = ForceFromPotential("Harmonic")
        integrator = EM(force, 1.0, 0.001)

        trajs=run_trajectories(integrator; params = params, init_conds_args = init_conds_args)
        list_files=filter(s -> occursin(r"trajectory_.*.dat",s),readdir())

        @test length(list_files) == 10
        for n in 1:length(list_files)
            data=readdlm(list_files[n])
            @test data[end,1] ≈ trajs[n].time[end]
            @test data[end,2:end] ≈ trajs[n].xt[end,:]
        end
        #=  Vérifier ici si on a bien sauvegardé ce qu'on voulait
        Donc on teste:
        - Qu'il y a bien 10 fichier différent
        - Que les fichiers font bien la taille attendu
        - Que l'état final du fichier est le même que la traj
        =#
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
