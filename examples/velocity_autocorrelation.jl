#=
main.jl

Compute velocity autocorrelation function for various set of integrators and compare it to expected value
=#

using Plots
using LangevinIntegrators
using StatsBase

let
    force=ForceFromPotential("Harmonic")

    params=TrajsParams(n_steps = 10^5, n_trajs = 10, n_save_iters = 1)
    plot()
    for int_class in [ABOBA,BAOAB,BBK,GJF]
        println(String(Symbol(int_class)))
        integrator=int_class(force, 1.0 , 1.0, 1.0, 1e-3, 1)

        trajs=run_trajectories(integrator; params = params, save_velocity=true)
        #Compute velocity auto correlation and plot it
        corr=zeros(div(params.n_steps,2)+1)
        for trj in trajs
            corr.+=autocov(trj.vt,0:div(params.n_steps,2))[:,1]/params.n_trajs
        end
        tps = LinRange(0,div(params.n_steps,2)*1e-3,div(params.n_steps,2)+1)
        plot!(tps, corr,label=String(Symbol(int_class)))
    end
    xlabel!("t")
    ylabel!("<v_0,v_t>")


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
    # tps = LinRange(0,10,10^4)
    # plot!(tps, traj[:,1],label="Traj")
    # xlabel!("t")
    # ylabel!("x")
end
