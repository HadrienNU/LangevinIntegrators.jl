#=
main.jl

Compute diffusion coefficients for various set of integrators and compare it to expected value
=#

using Plots
using LangevinIntegrators
using Distributions

let
    force=ForceFromPotential("Harmonic")
    γ = 1.0
    β = 1.0
    length_trj = 100.0
    t_range = LinRange(1e-3,1e-1,5)
    # inspectdr()
    plot(t_range, ones(length(t_range)),xaxis=:log, label="Theory")
    for int_class in [ABOBA,BAOAB,BBK,GJF]

        var_position=zeros(size(t_range))
        err_var=zeros(size(t_range)
        for (n,Δt) in enumerate(t_range)
            params=TrajsParams(n_steps = floor(length_trj/Δt), n_trajs = 250, n_save_iters = floor(0.01*length_trj/Δt))
            println(String(Symbol(int_class))," ", Δt)

            integrator=int_class(force, β , γ, 1.0, Δt, 1) # Change also initial conidition
            trajs=run_trajectories(integrator; params = params,to_save=["v"])
            x = vcat([trj.xt[1][:,1].^2 for trj in trajs]...)
            var_position[n] = mean(x)
            err_var[n] =quantile(TDist(length(x)-1), 1 - 0.05/2) * std(x)/sqrt(length(x))
        end
         plot!(t_range, var_position, y_err = err_var,xaxis=:log,marker=(:circle,1),label=String(Symbol(int_class)))
    end
    xlabel!("Δt")
    ylabel!("<x_t^2>")
end
