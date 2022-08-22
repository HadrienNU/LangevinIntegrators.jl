using LangevinIntegrators
using Test
# using Random

@testset "LangevinIntegrators.jl" begin
    # Write your tests here.
    @testset "Forces" begin
        # From potential
        force=ForceFromPotential("Harmonic")
        @test force_eval(force,[0.0])[1] ≈ 0.0
        @test force_eval(force,[1.0])[1] ≈ -1.0

        force=ForceFromPotential("DoubleWell")
        @test force_eval(force,[0.0])[1] ≈ 0.0
        @test force_eval(force,[1.0])[1] ≈ 0.0

        # force=ForceFromPotential("Muller")
        # @test force_eval(force,[0.0 0.0]) ≈ 0.0
        # @test force_eval(force,[1.0 1.0])[1] ≈ 0.0

        #Forces from ApproxFun
        force=ForceFromBasis("Taylor",[0.0 -2.0])
        @test force_eval(force,[0.0])[1] ≈ 0.0
        @test force_eval(force,[1.0])[1] ≈ -2.0
        #Forces from BSplinesKit
    end

    @testset "integrators_overdamped" begin
        force=ForceFromPotential("Harmonic")

        integrator=EM(force,1.0,1e-3)
        state = InitState!([0.0], integrator)
        UpdateState!(state, integrator)

    end


    @testset "integrators_inertial" begin
        force=ForceFromPotential("Harmonic")

        integrator=BBK(force, 1.0, 1.0, 1.0, 1e-3)
        state = InitState!([0.0],[0.0], integrator)
        UpdateState!(state, integrator)

        integrator=GJF(force, 1.0, 1.0, 1.0, 1e-3)
        state = InitState!([0.0],[0.0], integrator)
        UpdateState!(state, integrator)

        integrator=ABOBA(force, 1.0, 1.0, 1.0, 1e-3)
        state = InitState!([0.0],[0.0], integrator)
        UpdateState!(state, integrator)

        integrator=BAOAB(force, 1.0, 1.0, 1.0, 1e-3)
        state = InitState!([0.0],[0.0], integrator)
        UpdateState!(state, integrator)

        integrator=Verlet(force, 1.0, 1e-3)
        state = InitState!([0.0],[0.0], integrator)
        UpdateState!(state, integrator)

    end

    @testset "integrators_hidden" begin
        force=ForceFromPotential("Harmonic")
        params=LangevinParams()
        integrator=EM_Hidden(force,[[1.0,1.0] [-1.0,2.0]],[[1.0,0.0] [0.0,1.0]],1e-3,1)
        state = InitState!([0.0],[0.0],[0.0], integrator)
        UpdateState!(state, integrator)

        integrator=ABOBA_Hidden(force,[[1.0,1.0] [-1.0,2.0]],[[1.0,0.0] [0.0,1.0]],1e-3,1)
        state = InitState!([0.0],[0.0],[0.0], integrator)
        UpdateState!(state, integrator)

    end

end


# @testset "Plumed" begin
#     using LangevinIntegrators.Plumed
#     delta_t=2e-3
#     dim=4
#     temperature=1.0
#     energy=0.0
#     pos = [1.0,2.0,3.0,4.0] #Vector{Float}(undef,dim)
#     forces = [0.1,1.5,2.5,3.5] #Vector{Float}(undef,dim)
# end
#
# #Tout ce qui est en bas devrait aller dans les tests
# using .Plumed
#
#
# #Quelques infos à changer
# delta_t=2e-3
# dim=4
# temperature=1.0
#
# # A passer, en vrai il faudrait plutôt fournir pour l'energie une fonction callback qui ne calculera l'energie que si nécessaire.
# energy=0.0
# pos = [1.0,2.0,3.0,4.0] #Vector{Float}(undef,dim)
# forces = [0.1,1.5,2.5,3.5] #Vector{Float}(undef,dim)
#
# plumed_struct=plumed("../examples/plumed.dat","../examples/plumed.log",dim,delta_t,temperature)
#
# for n in 1:10
#     plumed_step(plumed_struct,n,pos,forces,energy)
#     println(forces)
# end
#
# plumed_finalize(plumed_struct)
