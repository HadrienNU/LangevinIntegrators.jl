using LangevinIntegrators.Plumed
# using LangevinIntegrators
# using Test

@testset "Plumed" begin

    # Testing first plumed by itself
    pos = [1.0, 2.0, 3.0, 4.0] #Vector{Float}(undef,dim)
    dim = length(pos)
    forces = zeros(dim)

    plumedFix = plumed("test_plumed.dat", "plumed.log", dim, 2e-3; temperature = 1.0)

    init_fix(plumedFix)

    for n = 1:10
        apply_fix!(plumedFix, pos, forces; step = n)
    end

    @test forces ≈ [-100.0, 0.0, 0.0, 0.0]

    @test plumedFix.bias_energy ≈ 10.0

    close_fix(plumedFix)

    @test isfile("plumed.log")
    @test isfile("colvar")

    # A new instance
    init_fix(plumedFix)
    @test isfile("bck.0.plumed.log")
    close_fix(plumedFix)


    # Testing in conjunction with force
    plumedFix = plumed("test_plumed_committor.dat", "plumed.log", 1, 0.001; temperature = 1.0)
    #Put it into a force and a integrator and run a trajectories

    init_conds_args = Dict("position" => Dict("type" => "Cste"))
    force = ForceFromPotential("Harmonic")
    addFix!(force, plumedFix)
    params = TrajsParams(; n_steps = 10^3)

    integrator = EM(force, 1.0, 0.001)
    init_conds = initialize_initcond(integrator, init_conds_args)

    state = InitState(integrator, init_conds)
    time, stop = run_trajectory!(state, integrator; params = params)

    # At the end, we should check that the colvar file contains the same value that the final state
    colvar = readdlm("colvar", comments = true)
    @test colvar[1, 1] ≈ integrator.Δt
    @test colvar[2, 1] ≈ 2*integrator.Δt # This check that plumed is not called twice in the same time step
    @test colvar[end, 1] ≈ time
    @test round(state.x[1], digits = 4) ≈ round(colvar[end, 2], digits = 4) # The precison  is lower in colvar file
    @test round(plumedFix.bias_energy, digits = 4) ≈ round(colvar[end, 3], digits = 4) # The precison  is lower in colvar file

    # We clean the test space by deleting the new file
    rm("colvar")
    rm("bck.0.colvar")
    rm("bck.1.colvar")
    rm("plumed.log")
    rm("bck.0.plumed.log")
    rm("bck.1.plumed.log")



    force = ForceFromPotential("Harmonic")
    integrator = EM(force, 1.0, 0.001)
    addPlumed!(integrator,"test_plumed_committor.dat", "plumed.log")
    @test length(integrator.force.fixes) == 1
end
