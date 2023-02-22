#=
main.jl

Compute diffusion coefficients for various set of integrators and compare it to expected value
=#

using Plots
using LangevinIntegrators
using Distributions

let
    force=ForceFromPotential("Flat")
    γ = 1.0
    β = 1.0
    length_trj = 100.0
    t_range = LinRange(1e-3,1e-1,5)
    # inspectdr()
    plot(t_range, ones(length(t_range)),xaxis=:log, label="Theory")
    for int_class in [ABOBA,BAOAB,BBK,GJF]

        diff_coeff=zeros(size(t_range))
        err_diff = zeros(size(t_range))
        for (n,Δt) in enumerate(t_range)
            params=TrajsParams(n_steps = floor(length_trj/Δt), n_trajs = 250, n_save_iters = floor(length_trj/Δt))
            println(String(Symbol(int_class))," ", Δt)

            integrator=int_class(force, β , γ, 1.0, Δt, 1)
            trajs=run_trajectories(integrator; params = params)
            end_points = [(trj.xt[1][end,1])^2/(2*params.n_steps*Δt) for trj in trajs]
            diff_coeff[n] = mean(end_points)
            err_diff[n] =quantile(TDist(params.n_trajs-1), 1 - 0.05/2) * std(end_points)/sqrt(params.n_trajs)
        end
         plot!(t_range, diff_coeff, y_err = err_diff,xaxis=:log,marker=(:circle,1),label=String(Symbol(int_class)))  # L'erreur est alors sqrt(var_diff)/params.n_trajs
    end
    xlabel!("Δt")
    ylabel!("<(x_t-x_0)^2>")
end
