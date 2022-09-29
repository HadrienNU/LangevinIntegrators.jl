#=
main.jl

Compute velocity autocorrelation function for various set of integrators and compare it to expected value
=#

using Plots
using LangevinIntegrators

let
    force=ForceFromPotential("Harmonic")

    params=TrajsParams()
    for int_class in [ABOBA,BAOAB,BBK,GJF]
        integrator=int_class(force, 1.0 , 1.0, 1.0, 1e-3, 1)
        state=InitState(integrator, params)
        traj=run_trajectories(integrator; params = params)
        #Compute velocity auto correlation and plot it
    end


    # for int_class in [BBK_Kernel,GJF_Kernel,EM_Kernel]
    #     integrator=int_class(force, 1.0 , 1.0, 1.0, 1e-3, 1)
    #     state=InitState(integrator, params)
    #     traj=run_trajectories(integrator; params = params)
    #     #Compute velocity auto correlation and plot it
    # end

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
