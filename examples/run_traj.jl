#=
main.jl

Read parameters file and launch the trajectories
=#

using Plots
include("../src/LangevinIntegrators.jl")
using .LangevinIntegrators


let
    # params,integrator=read_conf("onetraj.ini")
    integrator=EM(ForceFromPotential("Harmonic"),1.0,0.001)
    params=LangevinParams()
	state=init_trajectory(integrator; params = params)

	traj=run_trajectory!(state, integrator; params =params)
    println(state.x)


    # histogram([X[1] for X in traj],label="Traj",normalize=true);
    # qq = LinRange(-2,2,401)
    # plot!(qq, exp.(-β*V.(qq))/Z,label="Density")
    # xlabel!("x")
    # ylabel!("Frequency")
    tps = LinRange(0,10,10^4)
    plot!(tps, traj[:,1],label="Traj")
    xlabel!("t")
    ylabel!("x")
end
