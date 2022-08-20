using LangevinIntegrators
using Test

@testset "LangevinIntegrators.jl" begin
    # Write your tests here.
end


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
