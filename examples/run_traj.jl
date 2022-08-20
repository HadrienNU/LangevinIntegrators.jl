#=
main.jl

Read parameters file and launch the trajectories
=#

using Plots
using LangevinIntegrators

# LangevinIntegrators.TimerOutputs.enable_debug_timings(Main)  #Enable timer
# LangevinIntegrators.TimerOutputs.reset_timer!(LangevinIntegrators.timer) # Reset timer
let
    # params,integrator=read_conf("onetraj.ini")
    # force=ForceFromPotential("Harmonic")
    coeffs=zeros(Float64,(1,2))
    coeffs[1,:]=[0,-2]
    force=ForceFromBasis("Taylor",coeffs)
    integrator=EM(force,1.0,0.001)
    params=LangevinParams()
	state=init_trajectory(integrator; params = params)

	traj=run_trajectory!(state, integrator; params =params)

    # LangevinIntegrators.TimerOutputs.reset_timer!(LangevinIntegrators.timer) # Reset timer
    state=init_trajectory(integrator; params = params)
	traj=run_trajectory!(state, integrator; params =params)
    println(state.x)

    # println(LangevinIntegrators.timer)


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
