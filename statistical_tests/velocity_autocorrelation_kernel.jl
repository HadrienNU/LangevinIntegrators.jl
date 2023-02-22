#=
main.jl

Compute velocity autocorrelation function for various set of integrators and compare it to expected value
=#

using Plots
using LangevinIntegrators
using StatsBase

let
    force=ForceFromPotential("Harmonic")

    params=TrajsParams(n_steps = 5*10^5, n_trajs = 2, n_save_iters = 1)
    Δt = 1e-4
    tps = LinRange(0,div(params.n_steps,2)*Δt,div(params.n_steps,2)+1)
    kernel= exp.(-20.0*LinRange(0,500*1e-3, 500))
    β = 1.0
    γ=kernel[1]
    ω = sqrt(1-γ^2/4)
    # inspectdr()
    plot(tps,exp.(-0.5*γ*tps).*(cos.(ω*tps)-0.5*γ/ω*sin.(ω*tps))/β, label="Theory")
    for int_class in [BBK_Kernel,GJF_Kernel,EM_Kernel]
        println(String(Symbol(int_class)))
        integrator=int_class(force, β , kernel, 1.0, Δt, 1)

        trajs=run_trajectories(integrator; params = params, save_velocity=true)
        #Compute velocity auto correlation and plot it
        corr=zeros(div(params.n_steps,2)+1)
        for trj in trajs
            corr.+=autocov(trj.vt[10^3:end],0:div(params.n_steps,2))[:,1]/params.n_trajs
        end

        plot!(tps, corr,label=String(Symbol(int_class)))
    end
    xlabel!("t")
    ylabel!("<v_0,v_t>")

end
